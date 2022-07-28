# Flask test app
This is a test flask app. You can input an [SMILES string](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system) and will receive an image of the 2D structure of the corresponding molecule.

Furthermore, this app implements an authentication procces around an already existing user database. This database is initially created with the "add_database_entry.py" file and creates an user "guest" with the password "topsecret".

It is possible to change the password. Furthermore, inactive user of the app will be signed out after 60 seconds.

# Caveats
In a previous version, this app was "ready to be deployed" as a Docker container. In the current version, this is no longer the case due to added preexisting user database. The purpose of this repository simply shifted in the meantime.