from app import create_app

app = create_app()
app.app_context().push()
app.run(host="0.0.0.0")
