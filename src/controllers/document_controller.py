from http import HTTPStatus
import os
from flask import Blueprint, jsonify, request, current_app
from flask_jwt_extended import jwt_required
from werkzeug.utils import secure_filename
from src.services.doc_verification_service import doc_verification_service
base_dir = os.path.abspath(os.path.dirname(__file__))

UPLOAD_FOLDER = os.path.join(base_dir, "static/uploads") 
document = Blueprint("document_controller",__name__,url_prefix="/api/v1/document/validate")
ALLOWED_EXTENSIONS = set(['png', 'jpg'])

def remove_all_image():
    #Remove all the photos after processing them
    for image in os.listdir(current_app.config['UPLOAD_FOLDER']):
        os.remove(os.path.join(current_app.config['UPLOAD_FOLDER'], image))


@document.post("/passport")
@jwt_required()
def verify_passport():
    
    ##Upload passport image
    if "file" not in request.files:
        return jsonify({
            'error': 'media not provide'
            }), HTTPStatus.BAD_REQUEST
    file = request.files["file"]
    
    if file == '':
        return jsonify({
            'error': 'no file selected'
        }), HTTPStatus.BAD_REQUEST
    
    if file: 
        filename = secure_filename(file.filename)
        file.save(os.path.join(current_app.config['UPLOAD_FOLDER'], filename))

    #Get uploaded passport image for processing
    passport_image = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)
    
    doc_verification = doc_verification_service(passport_image)
    forgery_result =  doc_verification.check_image_forgery()
    
    if(forgery_result):#80 being the lowest prediction  
        
        mrz_result = doc_verification.check_MRZ()
    
        if(mrz_result):
            mrz_data = doc_verification.extract_info()
            if(doc_verification.validate_passport_details()):
                return jsonify({
                    "authentic": True,
                    "forged": False,
                    "validDetails": True,
                    "passportData":{
                        'Nationality': mrz_data['nationality'],
                        'Given Name': mrz_data['names'],
                        'Surname': mrz_data['surname'],
                        'Passport type': mrz_data['type'],
                        'Date of birth': mrz_data['date_of_birth'],
                        'ID Number': mrz_data['personal_number'],
                        'Gender': mrz_data['sex'],
                        'Expiration date': mrz_data['expiration_date']
                    }
                }),HTTPStatus.OK
            else:
                return jsonify({
                    'authentic': True,
                    'forged': False,
                    'message': "Passport details not valid"
                }),HTTPStatus.OK
        else:
            return jsonify({
                'authentic': mrz_result,
                'message': "Could not detect passport Machine Readable Zone"
            }),HTTPStatus.OK
    else:
        return jsonify({
            'authentic': forgery_result,
            'forged': True,
            'message': "Passport image might have been tempered with"
        }),HTTPStatus.OK

