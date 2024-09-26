FastAPI Application Deployment Guide
This guide explains the process of deploying the FastAPI application to an AWS EC2 instance using GitHub Actions for continuous deployment.

Overview
The deployment is automated using a GitHub Actions workflow that triggers on every push to the main branch. The workflow securely connects to an EC2 instance via SSH, syncs the application code, installs required dependencies, and starts the FastAPI application.

Key Components
GitHub Actions: Automates the deployment process triggered by code changes.
AWS EC2: The application is hosted on an EC2 instance running Ubuntu.
FastAPI: The core of the application, running as a backend API.
Uvicorn: ASGI server for running FastAPI.
rsync: Used to transfer files from GitHub to the EC2 instance over SSH.
Prerequisites
Before deploying, ensure you have the following:

An EC2 instance running Ubuntu, with SSH access.
SSH key pair configured and available as GitHub secrets.
Python 3.7+ installed on the EC2 instance.
Uvicorn installed as the ASGI server for the FastAPI app.
Deployment Process
The deployment is handled by a GitHub Actions workflow (deploy.yml). Below is a step-by-step explanation of what each part of the workflow does.

GitHub Actions Workflow (deploy.yml)
The workflow is triggered on every push to the main branch. It consists of the following steps:

Checkout the code: Uses actions/checkout@v4 to clone the repository into the workflow environment.

Install SSH Key: Uses the shimataro/ssh-key-action@v2 to install an SSH key from the GitHub secrets (EC2_SSH_KEY) to securely access the EC2 instance.

Add Known Hosts: Adds the EC2 host to the known_hosts file by scanning it with ssh-keyscan, ensuring the server's SSH key is trusted.

Deploy with rsync: Uses rsync to synchronize the project files from the GitHub repository to the EC2 instance, excluding certain files and directories such as .git and .github.

Execute remote SSH commands: Runs commands on the EC2 instance via appleboy/ssh-action. This includes:

Navigating to the project directory.
Installing Docker and Docker Compose on the EC2 instance.
Setting up a Python virtual environment, installing dependencies, applying database migrations with Alembic, and running the application using Docker Compose.


Secrets
To securely connect and deploy to your EC2 instance, the following secrets must be set in your GitHub repository:

EC2_SSH_KEY: The private SSH key used to authenticate access to the EC2 instance.
EC2_USER: The username used to connect to the EC2 instance (e.g., ubuntu).
EC2_HOST: The public IP address or DNS of the EC2 instance.
KNOWN_HOSTS: The public SSH host keys for the EC2 instance to avoid prompt confirmation during SSH connection.

Running the Application
After deployment, the application will be running and accessible via the public IP of the EC2 instance on port 80:

http://ec2-44-223-106-249.compute-1.amazonaws.com

To check the logs or status of the running application, you can SSH into the EC2 instance and view the Uvicorn logs:

ssh ubuntu@<ec2-44-223-106-249.compute-1.amazonaws.com>
tail -f nohup.out

Conclusion
This automated deployment pipeline ensures that every push to the main branch results in an updated FastAPI application running on your AWS EC2 instance. The use of GitHub Actions simplifies the process by automating SSH connections, code synchronization, and application startup.
