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

Checkout the Code: The actions/checkout@v4 action is used to fetch the latest version of the code from the repository.

Install SSH Key: The ssh-key-action installs the private SSH key from GitHub secrets to securely access the EC2 instance.

Add Known Hosts: This step ensures that the EC2 host is recognized, avoiding any Host verification failed errors during SSH connections.

Sync Files with EC2 Using rsync: The rsync command transfers the code from the GitHub repository to the EC2 instance. It excludes unnecessary files like .git, .github, and the SSH key itself.

Execute Remote Commands on EC2: This step runs a set of commands over SSH to configure and run the FastAPI app:

Update the package list.
Install python3-pip and any required Python packages.
Start the FastAPI application using Uvicorn with the following command:
bash
Copy code
nohup uvicorn main:app --host 0.0.0.0 --port 8000 &

Secrets
To securely connect and deploy to your EC2 instance, the following secrets must be set in your GitHub repository:

EC2_SSH_KEY: The private SSH key used to authenticate access to the EC2 instance.
EC2_USER: The username used to connect to the EC2 instance (e.g., ubuntu).
EC2_HOST: The public IP address or DNS of the EC2 instance.
KNOWN_HOSTS: The public SSH host keys for the EC2 instance to avoid prompt confirmation during SSH connection.
Running the Application
After deployment, the application will be running and accessible via the public IP of the EC2 instance on port 8000:

http://ec2-44-223-106-249.compute-1.amazonaws.com:8000

To check the logs or status of the running application, you can SSH into the EC2 instance and view the Uvicorn logs:

ssh ubuntu@<ec2-44-223-106-249.compute-1.amazonaws.com>
tail -f nohup.out

Conclusion
This automated deployment pipeline ensures that every push to the main branch results in an updated FastAPI application running on your AWS EC2 instance. The use of GitHub Actions simplifies the process by automating SSH connections, code synchronization, and application startup.
