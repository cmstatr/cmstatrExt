# average_curve_lm produces expected output

    Code
      print(res)
    Output
      
      Range: ` Strain ` in  [ 0,  0.1409409 ]
      
      Call:
      average_curve_lm(data = ., coupon_var = Coupon, model = Stress ~ 
          I(Strain) + I(Strain^2) + I(Strain^3) + 0, n_bins = 100)
      
      Coefficients:
        I(Strain)  I(Strain^2)  I(Strain^3)  
             1173        -8762        20481  
      

---

    Code
      summary(res)
    Output
      
      Range: ` Strain ` in  [0,  0.1409409 ]
      n_bins =  100 
      
      Call:
      average_curve_lm(data = ., coupon_var = Coupon, model = Stress ~ 
          I(Strain) + I(Strain^2) + I(Strain^3) + 0, n_bins = 100)
      
      Residuals:
          Min      1Q  Median      3Q     Max 
      -2.3080 -0.4003 -0.1726  0.3103  2.2058 
      
      Coefficients:
                   Estimate Std. Error t value Pr(>|t|)    
      I(Strain)    1173.285      4.662  251.69   <2e-16 ***
      I(Strain^2) -8761.916    102.493  -85.49   <2e-16 ***
      I(Strain^3) 20480.874    537.832   38.08   <2e-16 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      
      Residual standard error: 0.7546 on 397 degrees of freedom
      Multiple R-squared:  0.9997,	Adjusted R-squared:  0.9997 
      F-statistic: 3.985e+05 on 3 and 397 DF,  p-value: < 2.2e-16
      

# average_curve_optim produces expected output

    Code
      print(res_opt)
    Output
      
      Range: ` Strain ` in  [ 0,  0.1409409 ]
      
      Call:
      average_curve_optim(data = pa12_tension, coupon_var = Coupon, 
          x_var = Strain, y_var = Stress, fn = function(strain, par) {
              sum(par * c(strain, strain^2, strain^3))
          }, par = c(c1 = 1, c2 = 1, c3 = 1), n_bins = 100)
      
      Parameters:
             c1        c2        c3 
       1173.285 -8761.916 20480.874 

