from os import getenv
from fastapi import FastAPI
from src.molecules.router import router as molecule_router
from logging_config import setup_logging

setup_logging()


app = FastAPI()


@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}

app.include_router(molecule_router)



