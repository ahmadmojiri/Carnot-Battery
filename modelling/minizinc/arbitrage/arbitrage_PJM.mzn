% Use this editor as a MiniZinc scratch book
int: N;

set of int: tsamp = 1..N+1;    %number of sample points
set of int: tint = 1..N;     %number of intervals

float: dt;        %time difference between sample points (also, interval length) (h)
float: eta_in;    %efficiency of electricity to heat
float: eta_out;   %efficiency of heat to electricity
float: Q_max;     %maximum stored energy (MWh-th)
float: Q_min;     %minimum stored energy (MWh-th)
float: Q_init;    %initial stored energy (MWh-th)
float: Pin_max;   %maximum power to storage (MW-th)
float: Pin_min;   %minimum power to storage (MW-th)
float: Pout_max;  %maximum power to grid (MW-th)
float: Pout_min;  %minimum power to grid (MW-th)
float: Pda_frac;  %maximum power allowed into DA market
float: Prt_frac;   %maximum power allowed into RT market
float: SH;        %storage hours (h)
float: loss;
array[tint] of float: c_da; 
array[tint] of float: c_rt;

% === VARIABLES ===

array[tint] of var float: Pout_da;    %power out (MW-e) into the day ahead market
array[tint] of var float: Pout_rt;    %power out (MW-e) into the real time market
array[tint] of var float: Pin_da;     %power in (MW-e) from the day ahead market
array[tint] of var float: Pin_rt;     %power in (MW-e) from the real time market
array[tsamp] of var float: Q;      %Energy in storage (MWh-th)
var float: obj;

obj = sum(i in tint)(c_da[i] * (Pout_da[i] - Pin_da[i]) + c_rt[i] * (Pout_rt[i] - Pin_rt[i]) );

%=== CONSTRAINTS ===

constraint Q[1] = Q_init;

constraint forall(i in tint)( (Pout_da[i]) >= Pout_min*Pda_frac);
constraint forall(i in tint)( (Pout_da[i]) <= Pout_max*Pda_frac);

constraint forall(i in tint)( (Pin_da[i]) >= Pin_min*Pda_frac);
constraint forall(i in tint)( (Pin_da[i]) <= Pin_max*Pda_frac);

constraint forall(i in tint)( (Pout_rt[i]) >= Pout_min*Prt_frac);
constraint forall(i in tint)( (Pout_rt[i]) <= Pout_max*Prt_frac);

constraint forall(i in tint)( (Pin_rt[i]) >= Pin_min*Prt_frac);
constraint forall(i in tint)( (Pin_rt[i]) <= Pin_max*Prt_frac);

constraint forall(i in tint)( (Pout_da[i] + Pout_rt[i]) >= Pout_min);
constraint forall(i in tint)( (Pout_da[i] + Pout_rt[i]) <= Pout_max);

constraint forall(i in tint)( (Pin_da[i] + Pin_rt[i]) >= Pin_min);
constraint forall(i in tint)( (Pin_da[i] + Pin_rt[i]) <= Pin_max);

constraint forall(i in tsamp)(Q[i] >= Q_min);
constraint forall(i in tsamp)(Q[i] <= Q_max);

constraint forall(i in tint)(Q[i+1] = Q[i] + dt * ( (Pin_da[i]+Pin_rt[i]) * eta_in - (Pout_da[i]+Pout_rt[i]) / eta_out - loss*Q[i]/24));

solve maximize obj;

%output [show(obj)]

output  [show(obj*dt)] ++ [";"] ++ [show(Pin_da)] ++ [";"] ++ [show(Pin_rt)] ++ [";"] ++ [show(Pout_da)] ++ [";"] ++ [show(Pout_rt)] ++ [";"] ++ [show(Q)];

