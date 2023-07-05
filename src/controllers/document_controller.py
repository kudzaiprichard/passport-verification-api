from http import HTTPStatus
import os
from flask import Blueprint, jsonify, request, current_app
from flask_jwt_extended import jwt_required
from werkzeug.utils import secure_filename
#from src.__init__ import app
base_dir = os.path.abspath(os.path.dirname(__file__))

UPLOAD_FOLDER = os.path.join(base_dir, "static/uploads") 
document = Blueprint("document_controller",__name__,url_prefix="/api/v1/document/validate")
ALLOWED_EXTENSIONS = set(['png', 'jpg'])

@document.post("/id-passport")
#@jwt_required()
def id_passport():

    if "file" not in request.files:
        return jsonify({
            'error': 'media not provide'
            }), HTTPStatus.BAD_REQUEST
    
    files = request.files.getlist("file")
    
    if not files or len(files) < 3:
        return jsonify({
                    'error': 'They seem to be some missing files'
                }), HTTPStatus.BAD_REQUEST
    else:   
        for file in files:
            if file == '':
                return jsonify({
                    'error': 'no file selected'
                }), HTTPStatus.BAD_REQUEST
    
    if files: 
        for file in files:
            filename = secure_filename(file.filename)
            file.save(os.path.join(current_app.config['UPLOAD_FOLDER'], filename))
    
    
    #Remove all the photos after processing them
    for image in os.listdir(current_app.config['UPLOAD_FOLDER']):
        os.remove(os.path.join(current_app.config['UPLOAD_FOLDER'], image)) 
        
    return jsonify({
        'msg':'media upload successfully'
    })


@document.post("/other")
# @jwt_required()
def other():
    if "file" not in request.files:
        return jsonify({
            'error': 'media not provide'
            }), HTTPStatus.BAD_REQUEST
    
    files = request.files.getlist("file")
    
    if  len(files) > 1 :
        return jsonify({
                    'error': f'Expected a single image got {len(files)}'
                }), HTTPStatus.BAD_REQUEST
    else:   
        for file in files:
            if file == '':
                return jsonify({
                    'error': 'no file selected'
                }), HTTPStatus.BAD_REQUEST
    
    if files: 
        for file in files:
            filename = secure_filename(file.filename)
            file.save(os.path.join(current_app.config['UPLOAD_FOLDER'], filename))
    
    #Remove all the photos after processing them
    for image in os.listdir(current_app.config['UPLOAD_FOLDER']):
        os.remove(os.path.join(current_app.config['UPLOAD_FOLDER'], image))
        
    return jsonify({
        'msg':'media upload successfully'
    })
    return "Document::other"

