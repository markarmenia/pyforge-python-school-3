from os import getenv
from fastapi import FastAPI
from src.molecules.router import router as molecule_router
from logging_config import setup_logging
from src.cache import redis_client
import redis.asyncio as redis

setup_logging()


app = FastAPI()

@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}

@app.get("/health")
async def health_check():
    client = await redis_client()
    try:
        await client.ping()
        return {"status": "healthy"}
    except redis.RedisError:
        return {"status": "unhealthy"}, 503
    finally:
        await client.close()

app.include_router(molecule_router)



