from flask import Blueprint, request, jsonify
from http import HTTPStatus
import validators
from werkzeug.security import check_password_hash, generate_password_hash
from src.models.user import User, db
from flask_jwt_extended import get_jwt_identity, jwt_required, create_access_token,create_refresh_token

auth = Blueprint("auth_controller",__name__,url_prefix="/api/v1/auth")

@auth.post("/register")
def register():
    company = request.json['company']
    email = request.json['email']
    password = request.json['password']
    
    if len(password) < 6:
        return jsonify({
            'error': "Password is too short"
            }),HTTPStatus.BAD_REQUEST
    
    if len(company) < 3:
        return jsonify({
            'error': "Company name too short"
            }),HTTPStatus.BAD_REQUEST
    
    if not company.isalnum() or " " in company:
        return jsonify({
            'error': "Company should be alphanumeric with no spaces"
            }),HTTPStatus.BAD_REQUEST
    
    if not validators.email(email):
        return jsonify({
            'error': "Email is not valid"
            }),HTTPStatus.BAD_REQUEST
    
    if User.query.filter_by(email=email).first() is not None:
        return jsonify({'error': "Email is already in use"}),HTTPStatus.CONFLICT
    
    pwd_hash = generate_password_hash(password)
    user = User(company=company, password=pwd_hash, email=email)
    db.session.add(user)
    db.session.commit()
    
    return jsonify({
            'message': "User created",
            'user': {
                "company" : company,
                "email"   : email
            }
        }),HTTPStatus.CREATED

@auth.post("/authenticate")
def authenticate():
    email = request.json.get('email','')
    password = request.json.get('password','')
    
    user = User.query.filter_by(email=email).first()
    
    if user:
        is_pass_correct = check_password_hash(user.password, password)
        if is_pass_correct:
            refresh = create_refresh_token(identity=user.id)
            access = create_access_token(identity=user.id)
            
            return jsonify({
                'tokens': {
                    'accessToken': access,
                    'refreshToken': refresh
                }
            })
    
    return jsonify({'error' : 'Wrong credentials'}), HTTPStatus.UNAUTHORIZED


@auth.post("/refresh-token")
@jwt_required(refresh=True)
def refresh_token():
    identity = get_jwt_identity()
    access = create_access_token(identity=identity)
    return jsonify({
        'accessToken': access
    }), HTTPStatus.OK