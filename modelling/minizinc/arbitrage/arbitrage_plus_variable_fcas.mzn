int: N;

set of int: tsamp = 1..N+1;    %number of sample points
set of int: tint = 1..N;     %number of intervals

float: dt;      %time difference between sample points (also, interval length) (h)
float: eta_in;    %efficiency of electricity to heat
float: eta_out;   %efficiency of heat to electricity
float: self_discharge;  %standby losses as a percentage per day
float: Q_max;                              %maximum stored energy (MWh-th)
float: Q_min;                                %minimum stored energy (MWh-th)
float: Q_init;                             %initial stored energy (MWh-th)
float: Pin_max;                            %maximum power to storage (MW-e)
float: Pin_min;                              %minimum power to storage (MW-e)
float: Pout_max;                    %maximum power to grid (MW-e)
float: Pout_min;                         %minimum power to grid (MW-e)
float: SH;                      %storage hours (h)
array[tint] of float: rrp; 
array[tint] of float: fcas_rrp;
array[tint] of float: fcas_disp;

% === VARIABLES ===
%array[tint] of float: energy_frac;
%array[tint] of float: fcas_frac;
array[tint] of var float: Penergy;    %power out (MW-e)
array[tint] of var float: Pfcas;    %power out (MW-e)
array[tint] of var float: Pin;     %power in (MW-e)
array[tsamp] of var float: Q;      %Energy in storage (MWh-th)
var float: obj;

obj = dt*sum(i in tint)(rrp[i] * (Penergy[i] - Pin[i]) + fcas_rrp[i]* Pfcas[i]);

%=== CONSTRAINTS ===

constraint Q[1] = 0;

constraint forall(i in tint)(Penergy[i] >= Pout_min);
constraint forall(i in tint)(Pfcas[i] <= Pout_max - Penergy[i]);

constraint forall(i in tint)(Pin[i] >= Pin_min);
constraint forall(i in tint)(Pin[i] <= Pin_max);

constraint forall(i in tint)(Pfcas[i] <= fcas_disp[i]);
constraint forall(i in tint)(Pfcas[i] >= Pout_min);

constraint forall(i in tsamp)(Q[i] >= Q_min);
constraint forall(i in tsamp)(Q[i] <= Q_max);

constraint forall(i in tint)(Q[i+1] = Q[i] + dt * ( Pin[i]*eta_in - (Penergy[i])/eta_out) - Q[i]*self_discharge*dt/24 );

solve maximize obj;

output  [show(Pin)] ++ [";"] ++ [show(Penergy)] ++ [";"] ++ [show(Pfcas)] ++ [";"] ++ [show(Q)] ++ [";"] ++ [show(obj)];

