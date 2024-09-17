import redis.asyncio as redis
from os import getenv
from dotenv import load_dotenv

# Load the environment variables from .env
load_dotenv()

redis_host = getenv("REDIS_HOST", "localhost")
redis_port = int(getenv("REDIS_PORT", 6379))

print(f"Connecting to Redis at {redis_host}:{redis_port}")


async def redis_client():
    return redis.from_url(
        f"redis://{redis_host}:{redis_port}",
        decode_responses=True
    )