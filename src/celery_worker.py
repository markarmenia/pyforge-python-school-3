import os
from celery import Celery

redis_host = os.getenv('REDIS_HOST', 'localhost')
redis_port = os.getenv('REDIS_PORT', 6379)

celery_app = Celery(
    __name__,
    broker=f"redis://{redis_host}:{redis_port}/0",
    backend=f"redis://{redis_host}:{redis_port}/0"
)

celery_app.autodiscover_tasks(['src.tasks'])
celery_app.conf.update(task_track_started=True, task_time_limit=300)