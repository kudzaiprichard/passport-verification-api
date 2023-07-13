import os
from flask import current_app, jsonify
from passporteye import read_mrz
import numpy as np
from keras.models import load_model
from ela import convert_to_ela_image
model = load_model("trained_model.h5")

class doc_verification_service:
    def __init__(self, passport):
        self.passport = passport
    
    def prepare_image(self,fname):
        image_size = (128, 128)
        return (
            np.array(convert_to_ela_image(fname[0], 90).resize(image_size)).flatten()/255.0
        ) 


    def predict_result(self,fname):
        
        model = load_model("trained_model.h5")
        # outputs for the classification
        class_names = ["Forged", "Authentic"] 

        test_image = self.prepare_image(fname)
        test_image = test_image.reshape(-1, 128, 128, 3)

        y_pred = model.predict(test_image)
        y_pred_class = round(y_pred[0][0])

        prediction = class_names[y_pred_class]
        if y_pred <= 0.5:
            confidence = f"{(1-(y_pred[0][0])) * 100:0.2f}"
        else:
            confidence = f"{(y_pred[0][0]) * 100:0.2f}"
        return (prediction, confidence)

    def check_image_forgery(self):
        (prediction, confidence) = self.predict_result(self.passport)
        
        if(confidence > 80): #and prediction > 0.6
            return True
        
        return False
    
    def check_MRZ(self):
        mrz = read_mrz(self.passport)
        if(mrz == None):
            return False
        return True

    def validate_passport_details(self):
        mrz = read_mrz(self.passport)
        data = mrz.to_dict()
        
        if(data['valid_number'] == True, 
            data['valid_expiration_date'] == True,
            data['valid_composite'] == True,
            data['valid_personal_number'] == True,
            data['valid_date_of_birth'] == True,
            data['valid_expiration_date'] == True
        ):
            return True
        return False
    
    def extract_info(self):
        mrz = read_mrz(self.passport)
        return mrz.to_dict()