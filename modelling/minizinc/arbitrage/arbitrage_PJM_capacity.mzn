% Use this editor as a MiniZinc scratch book
int: N;
int: D=365;

set of int: tsamp = 1..N+1;    %number of sample points
set of int: tint = 1..N;     %number of intervals
set of int: tday = 1..D;

float: dt;      %time difference between sample points (also, interval length) (h)
float: eta_in;    %efficiency of electricity to heat
float: eta_out;   %efficiency of heat to electricity
float: Q_max;                              %maximum stored energy (MWh-th)
float: Q_min;                                %minimum stored energy (MWh-th)
float: Pin_max;                            %maximum power to storage (MW-th)
float: Pin_min;                              %minimum power to storage (MW-th)
float: SH;                               %storage hours (h)
float: loss;


array[tint] of float: c_rt;
array[tint] of float: Pout_rt;
array[tday] of int: Q_constraint;      % intervals indicating the start of day

% === VARIABLES ===

array[tint] of var float: Pin_rt;     %power in (MW-e) from the real time market
array[tsamp] of var float: Q;      %Energy in storage (MWh-th)
var float: obj;

obj = sum(i in tint)(c_rt[i] * (Pout_rt[i] - Pin_rt[i]));

%=== CONSTRAINTS ===

constraint forall(i in tint)( (Pin_rt[i]) >= Pin_min ); 
constraint forall(i in tint)( (Pin_rt[i]) <= Pin_max );


constraint forall(i in Q_constraint)( Q[i] >= 0.95*Q_max );


constraint forall(i in tsamp)(Q[i] >= Q_min);
constraint forall(i in tsamp)(Q[i] <= Q_max);

constraint forall(i in tint)(Q[i+1] = Q[i] + dt * ( (Pin_rt[i]) * eta_in - (Pout_rt[i]) / eta_out - loss*Q[i]/24));

solve maximize obj;

%output [show(Pout_rt)]

output  [show(obj*dt)] ++ [";"] ++ [show(Pin_rt)] ++ [";"] ++ [show(Pout_rt)] ++ [";"] ++ [show(Q)];

