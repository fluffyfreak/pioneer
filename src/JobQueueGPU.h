// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef GPUJOBQUEUE_H
#define GPUJOBQUEUE_H

#include <deque>
#include <vector>
#include <map>
#include <string>
#include <memory>
#include "SDL_thread.h"

class JobQueueGPU;
class JobRunnerGPU;
class JobHandleGPU;

// represents a single unit of work that you want done
// subclass and implement:
//
// OnRun: called from main thread, and does the actual stuff you want done.
//        store all your data in the object itself.
//
// OnFinish: called from the main thread once the GPU completes the job.
//           this is where you deliver the results from the GPU
//
// OnCancel: optional. called from the main thread to tell the job that its
//           results are not wanted. it should arrange for OnRun to return
//           as quickly as possible. OnFinish will not be called for the job
class JobGPU {
public:
	JobGPU() : cancelled(false), m_handle(nullptr) {}
	virtual ~JobGPU();

	JobGPU(const JobGPU&) = delete;
	JobGPU& operator=(const JobGPU&) = delete;

	virtual void OnRun() = 0;
	virtual void OnFinish() = 0;
	virtual void OnCancel() {}

private:
	friend class JobQueueGPU;
	friend class JobHandleGPU;
	friend class JobRunnerGPU;

	void UnlinkHandle();
	const JobHandleGPU* GetHandle() const { return m_handle; }
	void SetHandle(JobHandleGPU* handle) { m_handle = handle; }
	void ClearHandle() { m_handle = nullptr; }

	bool cancelled;
	JobHandleGPU* m_handle;
};


// a runner wraps a single thread, and calls into the queue when its ready for
// a new job. no user-servicable parts inside!
class JobRunnerGPU {
public:
	JobRunnerGPU(JobQueueGPU *jq);
	~JobRunnerGPU();
	void SetQueueDestroyed();
	void Process();

private:

	JobQueueGPU *m_jobQueue;

	JobGPU *m_job;

	bool m_queueDestroyed;
};

// This is the RAII handle for a queued JobGPU. A job is cancelled when the
// JobHandleGPU is destroyed. There is at most one JobHandleGPU for each JobGPU
// (non-queued GPUJobs have no handle). JobHandleGPU is not copyable only
// moveable.
class JobHandleGPU {
public:
	JobHandleGPU() : m_job(nullptr), m_queue(nullptr) { }
	JobHandleGPU(JobHandleGPU&& other);
	JobHandleGPU& operator=(JobHandleGPU&& other);
	~JobHandleGPU();

	JobHandleGPU(const JobHandleGPU&) = delete;
	JobHandleGPU& operator=(const JobHandleGPU&) = delete;

	bool HasGPUJob() const { return m_job != nullptr; }
	JobGPU* GetGPUJob() const { return m_job; }

private:
	friend class JobQueueGPU;
	friend class JobGPU;
	friend class JobRunnerGPU;

	JobHandleGPU(JobGPU* job, JobQueueGPU* queue);
	void Unlink();

	JobGPU* m_job;
	JobQueueGPU* m_queue;
};

// the queue management class. create one from the main thread, and feed your
// jobs do it. it will take care of the rest
class JobQueueGPU {
public:
	JobQueueGPU();
	~JobQueueGPU();

	// call from the main thread to add a job to the queue. the job should be
	// allocated with new. the queue will delete it once its its completed
	JobHandleGPU Queue(JobGPU *job);

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
	void Cancel(JobGPU *job);

	// call from the main loop. this will call OnFinish for any finished jobs,
	// and then delete all finished and cancelled jobs. returns the number of
	// finished jobs (not cancelled)
	Uint32 ProcessGPUJobs();

private:
	friend class JobRunnerGPU;
	JobGPU *GetGPUJob();
	void Finish(JobGPU *job);

	std::deque<JobGPU*> m_queue;

	std::deque<JobGPU*> m_finished;

	std::unique_ptr<JobRunnerGPU> m_runner;

	bool m_shutdown;
};

#endif
