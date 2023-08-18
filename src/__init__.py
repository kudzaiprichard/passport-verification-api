from flask import Flask
from src.controllers.auth_controller import auth
from src.controllers.document_controller import document
import os
from src.models.user import db
from flask_jwt_extended import JWTManager

base_dir = os.path.abspath(os.path.dirname(__file__))
def create_app(test_config=None):
    app = Flask(
        __name__,
        instance_relative_config=True
        )
    
    if test_config is None:
        app.config.from_mapping(
            SECRET_KEY = os.environ.get("SECRET_KEY"),
            SQLALCHEMY_DATABASE_URI=os.environ.get("SQLALCHEMY_DB_URI"),
            JWT_SECRET_KEY = os.environ.get('JWT_SECRET_KEY'),
            UPLOAD_FOLDER = os.path.join(base_dir, "static\\uploads") 
        )
        app.config.from_mapping(test_config)

    db.app = app
    db.init_app(app)
    JWTManager(app)
    app.register_blueprint(auth)
    app.register_blueprint(document)
    
    with app.app_context():
        db.drop_all()
        db.create_all() 
    return app