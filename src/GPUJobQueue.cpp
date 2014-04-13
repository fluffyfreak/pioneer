// Copyright © 2008-2014 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "GPUJobQueue.h"
#include "StringF.h"
#include "utils.h"

void GPUJob::UnlinkHandle()
{
	if (m_handle)
		m_handle->Unlink();
}

//virtual
GPUJob::~GPUJob()
{
	UnlinkHandle();
}

GPUJobRunner::GPUJobRunner(GPUJobQueue *jq) :
	m_jobQueue(jq),
	m_job(nullptr),
	m_queueDestroyed(false)
{
}

GPUJobRunner::~GPUJobRunner()
{
	// if we have a job running, cancel it. the worker will return it to the
	// finish queue, where it will be deleted later, so we don't need to do that
	if (m_job) {
		m_job->UnlinkHandle();
		m_job->OnCancel();
	}
}

void GPUJobRunner::Process()
{
	GPUJob *job;

	// Lock to prevent destruction of the queue while calling GetGPUJob.
	if (m_queueDestroyed) {
		return;
	}
	job = m_jobQueue->GetGPUJob();

	std::unique_ptr<MsgTimer> timer;
	if(job)
		timer.reset(new MsgTimer);

	Profiler::Timer jobTimes;
	jobTimes.Start();

	while (job) {
		// record the job so we can cancel it in case of premature shutdown
		m_job = job;

		// run the thing
		job->OnRun();

		// Lock to prevent destruction of the queue while calling Finish
		if (m_queueDestroyed) {
			return;
		}
		m_jobQueue->Finish(job);

		m_job = nullptr;

		// get a new job. this will block normally, or return null during
		// shutdown 
		if (m_queueDestroyed) {
			return;
		}

		jobTimes.SoftStop();
		const double ms = jobTimes.millicycles();
		if( ms > 1.0 / 60.0 ) {
			job = nullptr;
			break;
		} else {
			job = m_jobQueue->GetGPUJob();
		}
	}

	if(timer.get())
		timer->Mark("GPUJobRunner::Process()");
}

void GPUJobRunner::SetQueueDestroyed()
{
	m_queueDestroyed = true;
}


GPUJobHandle::GPUJobHandle(GPUJob* job, GPUJobQueue* queue) : m_job(job), m_queue(queue)
{
	assert(!m_job->GetHandle());
	m_job->SetHandle(this);
}

void GPUJobHandle::Unlink()
{
	if (m_job) {
		assert(m_job->GetHandle() == this);
		m_job->ClearHandle();
	}
	m_job = nullptr;
	m_queue = nullptr;
}

GPUJobHandle::GPUJobHandle(GPUJobHandle&& other) : m_job(other.m_job), m_queue(other.m_queue)
{
	if (m_job) {
		assert(m_job->GetHandle() == &other);
		m_job->SetHandle(this);
	}
	other.m_job = nullptr;
	other.m_queue = nullptr;
}

GPUJobHandle& GPUJobHandle::operator=(GPUJobHandle&& other)
{
	if (m_job && m_queue)
		m_queue->Cancel(m_job);
	m_job = other.m_job;
	m_queue = other.m_queue;
	if (m_job) {
		assert(m_job->GetHandle() == &other);
		m_job->SetHandle(this);
	}
	other.m_job = nullptr;
	other.m_queue = nullptr;
	return *this;
}

GPUJobHandle::~GPUJobHandle()
{
	if (m_job && m_queue) {
		m_queue->Cancel(m_job);
	} else {
		Unlink();
	}
}


GPUJobQueue::GPUJobQueue() :
	m_shutdown(false)
{
	m_runner.reset( new GPUJobRunner(this) );
}

GPUJobQueue::~GPUJobQueue()
{
	// flag shutdown. protected by the queue lock for convenience in GetGPUJob
	m_shutdown = true;

	// Flag each job runner that we're being destroyed (with lock so no one
	// else is running one of our functions). Both the flag and the mutex
	// must be owned by the runner, because we may not exist when it's
	// checked.
	m_runner->SetQueueDestroyed();

	// delete the runner.
	m_runner.reset();

	// delete any remaining jobs
	for (std::deque<GPUJob*>::iterator i = m_queue.begin(); i != m_queue.end(); ++i)
		delete (*i);
}

GPUJobHandle GPUJobQueue::Queue(GPUJob *job)
{
	GPUJobHandle handle(job, this);

	// push the job onto the queue
	m_queue.push_back(job);
	return handle;
}

// called by the runner to get a new job
GPUJob *GPUJobQueue::GetGPUJob()
{
	// loop until a new job is available
	GPUJob *job = nullptr;
	while (!job) {
		// we're shutting down, so just get out of here
		if (m_shutdown) {
			return nullptr;
		}

		if (m_queue.empty()) {
			// no jobs
			return nullptr;
		} else {
			// got one, pop it and return it
			job = m_queue.front();
			m_queue.pop_front();
		}

	}

	return job;
}

// called by the runner when a job completes
void GPUJobQueue::Finish(GPUJob *job)
{
	m_finished.push_back(job);
}

// call OnFinish methods for completed jobs, and clean up
Uint32 GPUJobQueue::ProcessGPUJobs()
{
	PROFILE_SCOPED()
	Uint32 finished = 0;

	m_runner->Process();

	if( m_finished.empty() ) {
		return finished;
	}
	GPUJob *job = m_finished.front();
	m_finished.pop_front();

	assert(job);

	// if its already been cancelled then its taken care of, so we just forget about it
	if(!job->cancelled) {
		job->UnlinkHandle();
		job->OnFinish();
		finished++;
	}

	delete job;

	return finished;
}

void GPUJobQueue::Cancel(GPUJob *job) {

	// check the waiting list. if its there then it hasn't run yet. just forget about it
	for (std::deque<GPUJob*>::iterator i = m_queue.begin(); i != m_queue.end(); ++i) {
		if (*i == job) {
			i = m_queue.erase(i);
			delete job;
		}
	}

	// check the finshed list. if its there then it can't be cancelled, because
	// its alread finished! we remove it because the caller is saying "I don't care"
	for (std::deque<GPUJob*>::iterator i = m_queue.begin(); i != m_queue.end(); ++i) {
		if (*i == job) {
			i = m_finished.erase(i);
			delete job;
			break;
		}
	}

	// its running, so we have to tell it to cancel
	job->cancelled = true;
	job->UnlinkHandle();
	job->OnCancel();
}
