import datetime
from flask_sqlalchemy import SQLAlchemy

db = SQLAlchemy()

class User(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    company = db.Column(db.String, unique=True, nullable=False)
    email = db.Column(db.String, nullable=False, unique=True)
    password = db.Column(db.Text(), nullable=False)
    created_at = db.Column(db.DateTime, default=datetime.datetime.now)
    
    
    def __repr__(self) -> str:
        return 'User>>>>{self.username}' 
