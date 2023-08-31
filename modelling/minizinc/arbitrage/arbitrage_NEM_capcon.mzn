% Use this editor as a MiniZinc scratch book
int: N;        %number of intervals in a year
int: S;        %number if intervals with c>strike price

set of int: tsamp = 1..N+1;    %number of sample points
set of int: tint = 1..N;     %number of intervals
set of int: tspike = 1..S;

float: dt;      %time difference between sample points (also, interval length) (h)
float: eta_in;    %efficiency of electricity to heat
float: eta_out;   %efficiency of heat to electricity
float: Q_init;                              %maximum stored energy (MWh-th)
float: Q_max;                              %maximum stored energy (MWh-th)
float: Q_min;                                %minimum stored energy (MWh-th)
float: Pin_max;                            %maximum power to storage (MW-th)
float: Pin_min;                              %minimum power to storage (MW-th)
float: Pout_max;                            %maximum power to storage (MW-th)
float: Pout_min;                              %minimum power to storage (MW-th)
float: SH;                               %storage hours (h)
float: loss;


array[tint] of float: c;
array[tspike] of int: pout_constraint;      % intervals indicating the start of day

% === VARIABLES ===

array[tint] of var float: Pin;     %power in (MW-e) from the real time market
array[tint] of var float: Pout;
array[tsamp] of var float: Q;      %Energy in storage (MWh-th)
var float: obj;

obj = sum(i in tint)(c[i] * (Pout[i] - Pin[i]));

%=== CONSTRAINTS ===

constraint forall(i in tint)( (Pin[i]) >= Pin_min ); 
constraint forall(i in tint)( (Pin[i]) <= Pin_max );


constraint forall(i in tint)( Pout[i] >= Pout_min ); 
constraint forall(i in tint)( Pout[i] <= Pout_max );

constraint forall(i in pout_constraint)( Pout[i] = Pout_max );

constraint Q[1] = Q_init*Q_max;
constraint forall(i in tsamp)(Q[i] >= Q_min);
constraint forall(i in tsamp)(Q[i] <= Q_max);

constraint forall(i in tint)(Q[i+1] = Q[i] + dt * ( (Pin[i]) * eta_in - (Pout[i]) / eta_out - loss*Q[i]/24));

solve maximize obj;

%output [show(Pout)]

output  [show(obj*dt)] ++ [";"] ++ [show(Pin)] ++ [";"] ++ [show(Pout)] ++ [";"] ++ [show(Q)] ++ [";"] ++ [show(pout_constraint)];

