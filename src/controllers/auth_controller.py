from flask import Blueprint

auth = Blueprint("auth_controller",__name__,url_prefix="/api/v1/auth")

@auth.post("/register")
def register():
    return "User Registration"

@auth.post("/authenticate")
def authenticate(email, password):
    return "User Authenticated"


@auth.post("/refresh-token")
def refresh_token(email, password):
    return "Refresh Token"