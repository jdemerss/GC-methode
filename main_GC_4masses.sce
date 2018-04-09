clc()

// Traitement du signal


data = fscanfMat('signal.txt')
acc_b = data(:,2)
t     = data(:,1)
 
fig1 = scf(1); 
plot(t,acc_b)



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%% 4DDL %%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


//%% Données SI
M = 1000;
H = 10;
E = 25000000000;
a = 0.35;   
b = 0.35;    

//%% Détermination des caractéristiques du bâtiment
I = b*a^3/12;
K = 3*E*I/H^3;

Kb1 = K;
Kb2 = K;
Ki = K*10000;
Ks = K*2;
Mb1 = M;
Mb2 = M;
Mi = M/100; 
Ms = M*10;


//%Newmark
dt = 0.01;
m = 10;
dt_1 = dt;
dt_2 = dt/m;
t_dt = 0:dt:t($);
t_dt_2 = 0:dt_2:t_dt($);
acc_b_dt = interp1(t,acc_b,t_dt_2,'linear');
Gamma = 1/2;
Beta = 1/4;
n = length(t_dt);


//=====================================================
//========= Résolution globale sol et structure =======
//=====================================================

A = [ acc_b_dt ; acc_b_dt; acc_b_dt; acc_b_dt/2];
K_ = [Kb1 -Kb1 0 0; -Kb1 Kb1+Kb2 -Kb2 0; 0 -Kb2 Kb2+Ki -Ki ; 0 0 -Ki (Ki+Ks)];
M_ = [Mb1 0 0 0; 0 Mb2 0 0; 0 0 Mi 0; 0 0 0 Ms];

nn = (n-1)*m+1;

[ Acc_N,V_N,U_N] = fct_Glob_DT_cst_4masses(nn,Gamma,Beta,dt_2,M_,K_,A)

//%%% Tracé des résultats

fig2 = scf(2)
scf(fig2)
subplot(221)
plot(t_dt_2',Acc_N(1,:)','b',t_dt_2',Acc_N(2,:)','r',t_dt_2', Acc_N(3,:)','g',t_dt_2', Acc_N(4,:)','k');
xtitle('Réponse en accélération','Time, t[s]','a [m/s²]')

subplot(222)
plot(t_dt_2',V_N(1,:)','b',t_dt_2',V_N(2,:)','r',t_dt_2', V_N(3,:)','g',t_dt_2', V_N(4,:)','k');
xtitle('Réponse en vitesse','Time, t[s]','v [m/s]')

subplot(223)
plot(t_dt_2',U_N(1,:)','b',t_dt_2',U_N(2,:)','r',t_dt_2', U_N(3,:)','g',t_dt_2', U_N(4,:)','k');
xtitle('Réponse en déplacement','Time, t[s]','u [m]')


//========================================================
//========== Problème de sous-structuration ==============
//========================================================

M_1 = [Mb1 0 0;0 Mb2 0;0 0 Mi/2];
M_2 = [Mi/2 0;0 Ms];
K_1 = [Kb1  -Kb1     0 ;..
      -Kb1 Kb1+Kb2 -Kb2;..
       0    -Kb2    Kb2];
K_2 = [Ki -Ki;-Ki Ki+Ks];

C_1 = [0 0 1];
C_2 = [-1 0];


//[ a_N,v_N,u_N] = fct_GC_DT_cst_4masses(n,Gamma,Beta,alpha,dt,M_1,M_2,K_1,K_2,A,C_1,C_2);
[ a_GC,v_GC,u_N_GC] = fct_GC_DT_nc_4masses(n,Gamma,Beta,dt_1,dt_2,M_1,M_2,K_1,K_2,A,C_1,C_2);
[ a_N,v_N,u_N] = fct_PH_DT_nc_4masses(n,Gamma,Beta,dt_1,dt_2,M_1,M_2,K_1,K_2,A,C_1,C_2,t_dt)



fig3 = scf(3)
scf(fig3)

subplot(221)
plot(t_dt_2',U_N(1,:)','b',t_dt_2',u_N_GC(1,:)','r',t_dt_2',u_N(1,:)','k');
xtitle('Réponse en déplacement 2eme etage','Time, t[s]','u [m]')

subplot(222)
plot(t_dt_2',U_N(2,:)','b',t_dt_2',u_N_GC(2,:)','r',t_dt_2',u_N(2,:)','k');
xtitle('Réponse en déplacement 1er etage','Time, t[s]','u [m]')

subplot(223)
plot(t_dt_2',U_N(3,:)','b',t_dt_2',u_N_GC(3,:)','r',t_dt_2',u_N(4,:)','k');
xtitle('Réponse en déplacement à l interface','Time, t[s]','u [m]')

subplot(224)
plot(t_dt_2',U_N(4,:)','b',t_dt_2',u_N_GC(5,:)','r',t_dt_2',u_N(5,:)','k');
xtitle('Réponse en déplacement du sol (m=10)','Time, t[s]','u [m]')

legends(['Résolution globale','GC','PH'],[2 5 1],opt="lr")


disp('fin')
