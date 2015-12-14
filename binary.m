# binary-mixture-adsorption-kinetics
# author: fliu000
# updated version

function binary_2(c1, c2, maxLoop)
    D   = 2.8 * 1e-6;
    x_m = 4.12 * 1e-10;  # saturation adsorption
    a1  = 2.9 * 1e-8;  # a1 is Langmuir von Szyskowski constant
    a2  = 2.97 * 1e-7;  # a2 is Langmuir von Szyskowski constant
    t   = 1;
    x   = (1 : maxLoop);
    y   = (1 : maxLoop);
    
    # temp var for estimations
    x_1 = (1 : maxLoop);
    y_1 = (1 : maxLoop);
    x_2 = (1 : maxLoop);
    y_2 = (1 : maxLoop);
    x_3 = (1 : maxLoop);
    y_3 = (1 : maxLoop);
    c1_s = (1 : maxLoop);
    c2_s = (1 : maxLoop);
    time_axis = (1 : maxLoop);
    
    maxIteration = 1000;
    
    x(1) = 0;
    y(1) = 0;
    
    c1_s(1) = a1 * x(1) / (x_m - x(1) - y(1));
    c2_s(1) = a2 * y(1) / (x_m - x(1) - y(1));
    
    x(2) = 2 * c1 * sqrt(D * t / pi);
    y(2) = 2 * c2 * sqrt(D * t / pi);
    
    c1_s(2) = a1 * x(2) / (x_m - x(2) - y(2));
    c2_s(2) = a2 * y(2) / (x_m - x(2) - y(2));
    
    time_axis(1) = 0;
    time_axis(2) = t;
    for n = 3 : maxLoop
        tn = (t * (n-1));
        time_axis(n) = tn;
        
        # Step 1: make the first estimate
        x_1(n) =  x(n-1) + (x(n-1) - x(n-2));
        y_1(n) =  y(n-1) + (y(n-1) - y(n-2));

        for iter = 1 : maxIteration
            # Step 2: compute the c_s(n)
            c1_s(n) = a1 * x_1(n) / (x_m - x_1(n) - y_1(n));
            c2_s(n) = a2 * y_1(n) / (x_m - x_1(n) - y_1(n));
            
            # Step 3: make the second estimate
            x_2(n) = 2 * sqrt(D / pi) * (c1 * sqrt(tn) - def_integral(c1_s, 1, floor(sqrt(tn)/t), t, n));
            y_2(n) = 2 * sqrt(D / pi) * (c2 * sqrt(tn) - def_integral(c2_s, 1, floor(sqrt(tn)/t), t, n));
            
            # Step 4: check whether second estimates are within 1%
            # of first estimate
            if (abs(x_2(n)/x_1(n) - 1) < 0.01 && abs(y_2(n)/y_1(n) - 1) < 0.01)
                # If yes, use second estimate as x(n), y(n)
                x(n) = x_2(n);
                y(n) = y_2(n);
                break;
            else 
                # If no, make the third estimate.
                x_3(n) = 1/2 * (x_1(n) + x_2(n));
                y_3(n) = 1/2 * (y_1(n) + y_2(n));
    
                # go to step 2.
                x_1(n) = x_3(n);
                y_1(n) = y_3(n);
            end
        end
    end
    
    clf;
    plot(time_axis, x);
    hold on
    plot(time_axis, y);
    hold on
    plot(time_axis, x+y, 'b-o');
end

function result = def_integral(val_array, startIndex, endIndex, delta, iter)
    # This computes the below
    # sum from i = startIndex to endIndex
    # val_array(iter - i) * delta / (2 * sqrt(delta * i))
    # which corresponses to the desired definite integral
    integral_result = 0;
    for i = startIndex : endIndex
        integral_result += val_array(iter - i) * delta / (2 * sqrt(delta * i));
    end
    result = integral_result;
end
