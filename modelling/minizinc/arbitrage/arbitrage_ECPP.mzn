int: N;
set of int: tsamp = 1..N+1;    %number of sample points
set of int: tint = 1..N;     %number of intervals

float: dt;         %time difference between sample points (also, interval length) (h)
float: R;          %ramp rate in MW/hr
float: eta_in;     %efficiency of electricity to heat
float: eta_out;    %efficiency of the power block
float: Q_max;      %maximum stored energy (MWh-th)
float: Q_min;      %minimum stored energy (MWh-th)
float: Q_init;     %initial stored energy (MWh-th)
float: Pin_max;    %maximum power to storage (MW-th)
float: Pin_min;    %minimum power to storage (MW-th)
float: Pout_max;   %maximum power to grid (MW-th)
float: Pout_min;   %minimum power to grid (MW-th)
float: SH;         %storage hours (h)
float: loss;
array[tint] of float: c; 

% === VARIABLES ===

array[tsamp] of var float: Pout;    %power out (MW-e)
array[tsamp] of var float: Pin;     %power in (MW-e)
array[tsamp] of var float: Q;      %Energy in storage (MWh-th)
var float: obj;

obj = sum(i in tint)(  dt*c[i] * (   0.5*(Pout[i+1]+Pout[i]) - 0.5*(Pin[i+1]+Pin[i])  )   );

%=== CONSTRAINTS ===
constraint Q[1] = Q_init;

constraint forall(i in tsamp)(Pout[i] >= Pout_min);
constraint forall(i in tsamp)(Pout[i] <= Pout_max);

constraint forall(i in tsamp)(Pin[i] >= Pin_min);
constraint forall(i in tsamp)(Pin[i] <= Pin_max);

constraint forall(i in tsamp)(Q[i] >= Q_min);
constraint forall(i in tsamp)(Q[i] <= Q_max);

constraint forall(i in tint)( Pout[i+1]<= (Pout[i]+dt*R)  );
constraint forall(i in tint)( Pout[i+1]>= (Pout[i]-dt*R)  );

constraint forall(i in tint)( Pin[i+1]<= (Pin[i]+dt*R)  );
constraint forall(i in tint)( Pin[i+1]>= (Pin[i]-dt*R)  );

                                                        

constraint forall(i in tint)(Q[i+1] = Q[i] + dt * (   0.5*(Pin[i]+Pin[i+1]) * eta_in -
                                                     Pout_max*( 0.962 + 2.085*( 0.5*(Pout[i]+Pout[i])/Pout_max-0.4)) - loss*Q[i]/24 )  );


solve maximize obj;

output  [show(Pin)] ++ [";"] ++ [show(Pout)] ++ [";"] ++ [show(Q)] ++ [";"] ++ [show(obj)] ;

