
//****************************************************************************
// Defining Variables
//****************************************************************************

var y c k n invest a w r log_y log_invest log_c log_n log_w log_y_n;
varexo e;

//****************************************************************************
// Defining Parameters
//****************************************************************************

parameters beta theta gamma alpha delta rho sigma;

//****************************************************************************
// Setting Parameter Values
//****************************************************************************

beta = 0.984;
theta = 3.48;
alpha = .667;
delta = .025;
rho = .979;
sigma = .01;
gamma = 1.004;

//****************************************************************************
// Model Equations
//****************************************************************************

model;
[name = 'Household problem FOC 1']
theta * c * 1 / (1-n) = w;

[name = 'Household problem FOC 2']
c(+1) / c = beta / gamma * (alpha* exp(a(+1)) *k(+1)^(alpha-1)*n(+1)^(1-alpha)+1-delta);

[name = 'Law of Motion of Capital'] 
gamma * k = (1-delta) * k(-1) + invest;

[name = 'Production Function']
y = exp(a) * k^(alpha) * (n)^(1-alpha);

[name = 'Resource Constraint']
y = c + invest;

[name = 'Firm FOC 1: Real Wage = Marginal Product of Labor']
w = (1-alpha) * exp(a) * k^alpha * n^(-alpha);

[name = 'Firm FOC 2: Annualized Real Interst Rate = Marginal Product of Capital']
r = 4 * alpha * exp(a) * k^(alpha-1) * n^(1-alpha);

[name = 'Exogenous TFP Process']
a = rho * a(-1) + e;

[name = 'Definition log output']
log_y = log(y);

[name = 'Definition log investment']
log_invest = log(invest);

[name = 'Definition log consumption']
log_c = log(c);

[name = 'Definition log hours']
log_n = log(n);

[name = 'Definition log wage']
log_w = log(w);

[name = 'Definition log output per hour worked']
log_y_n = log(y) - log(n);

end;


//****************************************************************************
// Computing the Steady State Given the Starting Values
//****************************************************************************

initval;
y=30;
k=400;
c=15;
n=0.5;
invest=10;
r=0;
w=10;
a=1;
e=0;
end;

resid(1);
steady;

//****************************************************************************
// Setting Shock Variance
//****************************************************************************

shocks;
var e = sigma^2;
end;

//****************************************************************************
// Checking Blanchard-Kahn Conditions
//****************************************************************************

check;

//****************************************************************************
// Computing
//****************************************************************************

stoch_simul(irf=40, hp_filter=1600) log_y log_c log_invest log_n log_y_n log_w r a;
