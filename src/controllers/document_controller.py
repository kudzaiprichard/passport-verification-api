from flask import Blueprint
from flask_jwt_extended import jwt_required

document = Blueprint("document_controller",__name__,url_prefix="/api/v1/document/validate")

@document.post("/id-passport")
@jwt_required()
def id_passport():
    return "Document::id-passport"

@document.post("/other")
@jwt_required()
def other():
    return "Document::other"

