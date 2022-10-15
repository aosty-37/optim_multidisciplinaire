%*******************************************************************************
%Cette fonction algo_SQP_BFGS implemente l'algorithme SQP (sans derivee)
%                                                                              *
%                                                                              *
%                                                                              *
%*******************************************************************************
%                                                    %**************************
%                                                    % PARAMETRES EN ENTREE    *
%                                                    %**************************
%
%            une_f                % une fonction dont on cherche un minimum    %
%
%            des_c                % les contraintes des probleme               %
%
%            un_x0                % un incontournable point initial  
%
%           un_lambda0            % une valeur intiale pour les
%           multiplicateurs de Lagrange.
%            
%            H0                   % une initialisation pour la matrice H_k     %
%
%            une_tol_h            % une tol. pertubation pour approche les
%            derivees (via diffenreces finies)
%
%            un_nit_max           % nombre maximum d'iterations autorisees     %
%                                 % risees                                     %
%            une_tol_x            % seuil de stationnarite des x_k             %
%            une_tol_f            % seuil de stationnarite des f_k             %
%            une_tol_g            % seuil validant EULER    en x_k             %
%                                   le sous-probl�me:
%                                   (1 exact, 2 pas de Cauchy, 3 CG tronqu�)   %
%            varargin             % la fonction du noyau MATLAB qui collecte   %
%                                 % les parametres additionnels necessaires au %
%                                 % calcul de une_f, donc un_gf eventuellement %
%                                 % donc une_hf eventuellement                 %
%
%                                                    %**************************
%                                                    % PARAMETRES EN SORTIE    *
%                                                    %**************************
%
%            x_opt                % la solution proposee par trust_region      %
%            f_opt                % une_f (x_opt, varargin{:})                 %
%            g_opt                % un_gf (x_opt, varargin{:})                 %
%            fin                  % la cause de l'arret de l'algorithme        %
%            nit                  % le nombre iterations de trust_region       %
%            f_count              % le nombre d'evaluations de une_f           %
%            g_count              % le nombre d'evaluations de un_gf
%
% Responsable: Y. Diouane (youssef.diouane@isae.fr) -- 
% (C) Institut Sup�rieur de l'A�ronautique et de l'Espace (ISAE-SUPAERO)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [                                                             ...
    x_opt  ,                                                           ...
    f_opt  ,                                                           ...
    fin    ,                                                           ...
    nit                                                                ...
    ]                                                                  ...
    = algo_SQP_BFGS(                                ...
    une_f          ,                                ...
    des_c          ,                                ...
    un_x0          ,                                ...
    un_lambda0     ,                                ...
    H0             ,                                ...
    une_tol_h      ,                                ...
    un_nit_max     ,                                ...
    une_tol_x      ,                                ...
    une_tol_f      ,                                ...
    une_tol_g                                       ...
    )

%**********************
%    CODE             *
%**********************

global f_count  ;                  % nombre     d'evaluations de la
% fonction a minimiser, sans pitie

global g_count  ;                  % nombre     d'evaluations du
% gradient g_f        , sans pitie
global h_count  ;                  % nombre     d'evaluations de la
% Hessienne h_f       , sans pitie
global c_count  ;                  % nombre     d'evaluations des
% constraintes a minimiser, sans pitie

global jc_count  ;                  % nombre     d'evaluations du
% gradient jac_c        , sans pitie
global hc_count  ;                  % nombre     d'evaluations de la
% Hessienne h_c       , sans pitie

f_count    = 0                                                             ;
g_count    = 0                                                             ;
h_count    = 0                                                             ;

c_count    = 0                                                             ;
jc_count   = 0                                                             ;
hc_count   = 0                                                             ;

% tempx sera desormais    %
% vecteur colonne         %
tempx           = un_x0'                                                   ;
templambda      = un_lambda0'                                              ;
tempx_pred      = un_x0';
templambda_pred = un_lambda0';

n               = length(tempx   )                                         ;
p               = length(templambda   )                                    ;

% Evaluer notre Lagrangien
fdex            = feval(une_f,tempx)                                       ;
cdex            = feval(des_c,tempx)                                       ;

% Evaluer le gradient de notre Lagrangian
gfdex = zeros(n, 1) ;
jcdex = zeros(n, 1) ;
E = zeros(n,n); % matrice de base
for k = 1: n
    E(k, k) = 1.;
end
for i = 1:n
    gfdex(i) = fdex - feval(une_f, tempx +une_tol_h*E(i,:));  
end
for j = 1:n
    jcdex(j) = cdex - feval(des_c, tempx +une_tol_h*E(j, :));  
end
gfdex = gfdex/une_tol_h;
jcdex = jcdex'/une_tol_h;

% Evaluer le Hessien de notre Lagrangian
H               = H0                                                       ;
k               = 0                                                        ;
fin             = 0                                                        ;
while(fin==0 && k < un_nit_max)
    % Evaluer notre Lagrangien et gradient du Lagrangien � xk_pred
    gfdex_pred          = gfdex                                ;
    jcdex_pred          = jcdex                                 ;                                                      ; 

    fdex            = feval(une_f,tempx)                                       ;
    cdex            = feval(des_c,tempx)                                       ;
    
    
    % Calcul des gradients avec des diff�rences finies
    for i = 1:n
        gfdex(i) = fdex - feval(une_f, tempx +une_tol_h*E(i,:));  
    end
    for j = 1:n
        jcdex(j) = cdex - feval(des_c, tempx +une_tol_h*E(j,:));  
    end
    gfdex = gfdex /une_tol_h;
    jcdex = jcdex /une_tol_h;

    % Calcul du Hessien du Lagrangien avec le BFGS
    sk_pred = tempx - tempx_pred ;
    yk_pred = (gfdex + templambda' * jcdex ) - (gfdex_pred +templambda_pred'* jcdex_pred);

    if(yk_pred'*sk_pred >0)
        H = H + ((yk_pred*yk_pred') .*  (1/ (yk_pred'*sk_pred))) - ((H* (sk_pred* sk_pred') * H) .*  (1/ (sk_pred'*H* sk_pred))) ;
    end
    
    % R�solution du probl�me quadratique
    LSk = [H jcdex';jcdex 0]        ;
    yk = [ gfdex ; cdex] ;
    res=linsolve(LSk, -yk);
    dk = res(1:n);
    templambda_pred = templambda ;
    templambda = res(n+1:n+p);
    
    if (norm(dk) <une_tol_x)
        fin = 4;
    end
    if (norm(gfdex) < une_tol_g)
        fin = 1;
    end
    if ((fdex - feval(une_f, tempx+dk)) < une_tol_f)
        fin = 2;
    end
    if(k==un_nit_max)
        fin =  3                                                           ;
    end
    tempx_pred = tempx;
    tempx = tempx +dk ;
    k = k +1 ;
end
x_opt    =   tempx                                                         ;
f_opt    =    fdex                                                         ;
nit      =   k                                                             ;
