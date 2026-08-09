"""Small thread-safe runtime caches used by the Flask app."""

import copy
import threading
import time
import uuid
from collections import OrderedDict


class _DesignResultCache:
    """LRU + TTL cache with the dict-like API used by design results."""

    def __init__(self, max_size: int, ttl_seconds: int):
        self._max = max_size
        self._ttl = ttl_seconds
        self._data: "OrderedDict[str, tuple]" = OrderedDict()
        self._lock = threading.Lock()

    def _now(self) -> float:
        return time.monotonic()

    def __setitem__(self, key, value):
        with self._lock:
            self._data.pop(key, None)
            self._data[key] = (self._now() + self._ttl, value)
            self._data.move_to_end(key)
            while len(self._data) > self._max:
                self._data.popitem(last=False)

    def get(self, key, default=None):
        with self._lock:
            entry = self._data.get(key)
            if entry is None:
                return default
            expire_at, value = entry
            if expire_at < self._now():
                del self._data[key]
                return default
            self._data.move_to_end(key)
            return value

    def __contains__(self, key):
        return self.get(key) is not None

    def __len__(self):
        with self._lock:
            return len(self._data)

    def clear(self):
        with self._lock:
            self._data.clear()

    def stats(self):
        return {"size": len(self), "max": self._max, "ttl": self._ttl}


class _ActivityValidationCache:
    """Thread-safe in-memory task state for the activity validation workbench."""

    def __init__(self, max_size: int, ttl_seconds: int):
        self._max = max_size
        self._ttl = ttl_seconds
        self._data: "OrderedDict[str, tuple]" = OrderedDict()
        self._lock = threading.Lock()

    def _now(self) -> float:
        return time.monotonic()

    def create(self, *, job_id: str, activities: list) -> str:
        task_id = uuid.uuid4().hex[:12]
        now = self._now()
        with self._lock:
            self._cleanup_locked(now)
            self._data[task_id] = (
                now + self._ttl,
                {
                    "task_id": task_id,
                    "job_id": job_id,
                    "activities": list(activities),
                    "status": "queued",
                    "stage": "queued",
                    "progress": 0.0,
                    "events": [],
                    "partial_results": [],
                    "result": None,
                    "artifacts": [],
                    "error": None,
                    "created_at": time.time(),
                    "updated_at": time.time(),
                },
            )
            self._data.move_to_end(task_id)
            while len(self._data) > self._max:
                self._data.popitem(last=False)
        self.event(task_id, "queued", "validation task queued", progress=0)
        return task_id

    def update(self, task_id: str, **fields) -> None:
        with self._lock:
            entry = self._data.get(task_id)
            if entry is None:
                return
            expire_at, payload = entry
            payload.update(fields)
            payload["updated_at"] = time.time()
            self._data[task_id] = (expire_at, payload)
            self._data.move_to_end(task_id)

    def event(self, task_id: str, stage: str, message: str, *, progress=None, **detail) -> None:
        with self._lock:
            entry = self._data.get(task_id)
            if entry is None:
                return
            expire_at, payload = entry
            event = {
                "time": time.time(),
                "stage": stage,
                "message": message,
            }
            if detail:
                event["detail"] = detail
            payload.setdefault("events", []).append(event)
            payload["events"] = payload["events"][-120:]
            payload["stage"] = stage
            if progress is not None:
                payload["progress"] = max(0.0, min(100.0, float(progress)))
            payload["updated_at"] = time.time()
            self._data[task_id] = (expire_at, payload)
            self._data.move_to_end(task_id)

    def get(self, task_id: str):
        now = self._now()
        with self._lock:
            entry = self._data.get(task_id)
            if entry is None:
                return None
            expire_at, payload = entry
            if expire_at < now:
                del self._data[task_id]
                return None
            self._data.move_to_end(task_id)
            return copy.deepcopy(payload)

    def _cleanup_locked(self, now: float) -> None:
        for key in list(self._data.keys()):
            expire_at, _payload = self._data[key]
            if expire_at < now:
                del self._data[key]
