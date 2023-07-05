from flask import Blueprint

document = Blueprint("document_controller",__name__,url_prefix="/api/v1/document/validate")

@document.post("/id-passport")
def id_passport(request):
    return "Document::id-passport"

@document.post("/other")
def other(request):
    return "Document::other"

