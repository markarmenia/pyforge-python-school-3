import redis.asyncio as redis
from os import getenv

redis_host = getenv("REDIS_HOST", "localhost")
redis_port = int(getenv("REDIS_PORT", 6379))

async def redis_client():
    return redis.from_url(
        f"redis://{redis_host}:{redis_port}",
        decode_responses=True
    )
