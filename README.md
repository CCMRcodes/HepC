Copyright (C) 2019 Waljee Lab, The Regents of the University of Michigan

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see the full license text file.


Overview

This repository contains code necessary to perform the data analysis described in the paper “Machine Learning Models to Predict Disease Progression Among Veterans with Hepatitis C Virus”. This paper developed and compared two machine learning algorithms, Cox model and boosted-survival tree model, to predict cirrhosis development in a large CHC-infected cohort using longitudinal data.


Usage

To run the code, the user needs to preprocess the data into the format according to information in the relevant R files.

•Model_boosting_cross-sectional.R: This code performs cross-sectional boosting survival model with 30 random splits of training/test data.


•model_boosting_longitudinal.R: This code performs longitudinal boosting survival model with 30 random splits of training/test data.


•model_cox_cross-sectional_longitudinal.R: This code performs cross-sectional and longitudinal Cox models with 30 random splits of training/test data.


•model_misclassification_table.R: This code computes for the misclassification table under one representative split. The split is selected to have the closest concordance to average in the boosting longitudinal model.


•plot_partial_dependence_boosting_cross-sectional.R: This code produces partial dependence plots for boosting cross-sectional model.


•plot_partial_dependence_boosting_longitudinal.R: This code produces partial dependence plots for boosting longitudinal model.


•plot_variable_importance_boosting.R: This code produces variable importance plots for boosting cross-sectional and longitudinal models.


•variable_selection_cox_cross-sectional.R: This code performs variable selection for cross-sectional Cox model with 30 random splits of training/test data.


•variable_selection_cox_longitudinal.R: This code performs variable selection for longitudinal Cox model with 30 random splits of training/test data.



Reference

Monica A. Konerman, Lauren A. Beste, Tony Van, Boang Liu, Xuefei Zhang, Ji Zhu, Sameer D. Saini, Grace L. Su, Brahmajee K. Nallamothu, George N. Ioannou, Akbar K. Waljee. “Machine Learning Models to Predict Disease Progression Among Veterans with Hepatitis C Virus.”
