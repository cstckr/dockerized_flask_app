# Flask test app
This is a test flask app ready to be deployed as an Azure web app (via an docker image). You can input an [SMILES string](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system) and will receive an image of the 2D structure of the corresponding molecule.

This app implements an authentication procces around an Azure SQL user database. Inactive users of the app will be signed out after 120 seconds.


# Deploy to Azure

- Create database in Azure
 - In Azure, go to the "SQL databases" service and create a SQL database
 - In this repository, create a folder "credentials", that contains a "database.py" file that holds corresponding string variables for "server", "database", "username", and "password". This is ommited in this repository for obvious reasons.
 - Run the file "manage_azure_database.py" to create a user table in the Azure database with an user "guest3" with password "123456"

- Create docker image and push in to the Azure container registry
 - Start Docker Desktop
 - Install the [Azure CLI](https://docs.microsoft.com/en-us/cli/azure/install-azure-cli)
 - In Azure, go to the "Container registries" service and create a registry
 - In a command terminal, go to your repository folder and run following commands:
        docker run -dit ubuntu sh
		docker image build -t docker-app .
		az login
		az acr login --name YOURLOGINSERVER.azurecr.io
		docker tag docker-app YOURLOGINSERVER.azurecr.io/testapp
		docker push YOURLOGINSERVER.azurecr.io/testapp
		
 - In your Azure container registry, make sure to enable "Admin user" in the "Update" tab 

- In Azure, go to the "App Services" service and create the app service by selecting your docker image from the Azure container registry.