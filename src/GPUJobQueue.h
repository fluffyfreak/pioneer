// Copyright © 2008-2014 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef GPUJOBQUEUE_H
#define GPUJOBQUEUE_H

#include <deque>
#include <vector>
#include <map>
#include <string>
#include <memory>
#include "SDL_thread.h"

class GPUJobQueue;
class GPUJobRunner;
class GPUJobHandle;

// represents a single unit of work that you want done
// subclass and implement:
//
// OnRun: called from worker thread, and does the actual stuff you want done.
//        store all your data in the object itself.
//
// OnFinish: called from the main thread once the worker completes the job.
//           this is where you deliver the results from the worker
//
// OnCancel: optional. called from the main thread to tell the job that its
//           results are not wanted. it should arrange for OnRun to return
//           as quickly as possible. OnFinish will not be called for the job
class GPUJob {
public:
	GPUJob() : cancelled(false), m_handle(nullptr) {}
	virtual ~GPUJob();

	GPUJob(const GPUJob&) = delete;
	GPUJob& operator=(const GPUJob&) = delete;

	virtual void OnRun() = 0;
	virtual void OnFinish() = 0;
	virtual void OnCancel() {}

private:
	friend class GPUJobQueue;
	friend class GPUJobHandle;
	friend class GPUJobRunner;

	void UnlinkHandle();
	const GPUJobHandle* GetHandle() const { return m_handle; }
	void SetHandle(GPUJobHandle* handle) { m_handle = handle; }
	void ClearHandle() { m_handle = nullptr; }

	bool cancelled;
	GPUJobHandle* m_handle;
};


// a runner wraps a single thread, and calls into the queue when its ready for
// a new job. no user-servicable parts inside!
class GPUJobRunner {
public:
	GPUJobRunner(GPUJobQueue *jq);
	~GPUJobRunner();
	void SetQueueDestroyed();
	void Process();

private:

	GPUJobQueue *m_jobQueue;

	GPUJob *m_job;

	bool m_queueDestroyed;
};

// This is the RAII handle for a queued GPUJob. A job is cancelled when the
// GPUJobHandle is destroyed. There is at most one GPUJobHandle for each GPUJob
// (non-queued GPUJobs have no handle). GPUJobHandle is not copyable only
// moveable.
class GPUJobHandle {
public:
	GPUJobHandle() : m_job(nullptr), m_queue(nullptr) { }
	GPUJobHandle(GPUJobHandle&& other);
	GPUJobHandle& operator=(GPUJobHandle&& other);
	~GPUJobHandle();

	GPUJobHandle(const GPUJobHandle&) = delete;
	GPUJobHandle& operator=(const GPUJobHandle&) = delete;

	bool HasGPUJob() const { return m_job != nullptr; }
	GPUJob* GetGPUJob() const { return m_job; }

private:
	friend class GPUJobQueue;
	friend class GPUJob;
	friend class GPUJobRunner;

	GPUJobHandle(GPUJob* job, GPUJobQueue* queue);
	void Unlink();

	GPUJob* m_job;
	GPUJobQueue* m_queue;
};

// the queue management class. create one from the main thread, and feed your
// jobs do it. it will take care of the rest
class GPUJobQueue {
public:
	GPUJobQueue();
	~GPUJobQueue();

	// call from the main thread to add a job to the queue. the job should be
	// allocated with new. the queue will delete it once its its completed
	GPUJobHandle Queue(GPUJob *job);

	// call from the main thread to cancel a job. one of three things will happen
	//
	// - the job hasn't run yet. it will never be run, and neither OnFinished nor
	//   OnCancel will be called. the job will be deleted on the next call to
	//   FinishGPUJobs
	//
	// - the job has finished. neither onFinished not onCancel will be called.
	//   the job will be deleted on the next call to FinishGPUJobs
	//
	// - the job is running. OnCancel will be called
	void Cancel(GPUJob *job);

	// call from the main loop. this will call OnFinish for any finished jobs,
	// and then delete all finished and cancelled jobs. returns the number of
	// finished jobs (not cancelled)
	Uint32 ProcessGPUJobs();

private:
	friend class GPUJobRunner;
	GPUJob *GetGPUJob();
	void Finish(GPUJob *job);

	std::deque<GPUJob*> m_queue;

	std::deque<GPUJob*> m_finished;

	std::unique_ptr<GPUJobRunner> m_runner;

	bool m_shutdown;
};

#endif
