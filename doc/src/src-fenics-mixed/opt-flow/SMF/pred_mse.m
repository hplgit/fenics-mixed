function mse = pred_mse(x,dmodel)

[f,df,mse] = predictor(x, dmodel);
% minimize - mse to maximize mse
mse = -mse;
