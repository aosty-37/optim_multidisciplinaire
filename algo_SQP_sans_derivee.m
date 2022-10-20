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
%                                   le sous-problème:
%                                   (1 exact, 2 pas de Cauchy, 3 CG tronqué)   %
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
% (C) Institut Supérieur de l'Aéronautique et de l'Espace (ISAE-SUPAERO)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [                                                             ...
    x_opt  ,                                                           ...
    f_opt  ,                                                           ...
    fin    ,                                                           ...
    nit                                                                ...
    ]                                                                  ...
    = algo_SQP_sans_derivee(                                ...
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

n               = length(tempx)                                         ;
p               = length(templambda)                                    ;

% Evaluer notre Lagrangien
fdex            = feval(une_f,tempx)                                       ;
cdex            = feval(des_c,tempx)                                       ;

% Evaluer le gradient de notre Lagrangian
gfdex = zeros(n, 1) ;
jcdex = zeros(p, n) ;
E = eye(n);

for i =1:n
    gfdex(i, 1) =(feval(une_f, tempx +une_tol_h*E(:,i)) - fdex)/ une_tol_h;    
    jcdex(:, i) =(feval(des_c, tempx +une_tol_h*E(:,i)) - cdex) /une_tol_h;
end

% Evaluer le Hessien de notre Lagrangian
H               = H0                                                       ;
k               = 0                                                        ;
fin             = 0                                                        ;

% Initialisation de sk_pred et yk_pred 
sk_pred = 0;
yk_pred = 0;
while(fin==0 && k < un_nit_max)
 
    if(yk_pred' * sk_pred >0)
        alpha = (1/ (yk_pred'*sk_pred))* (yk_pred*yk_pred') ; 
        beta = (1/ (sk_pred'*H* sk_pred)) *   (H* (sk_pred* sk_pred') * H);
        H = H  + alpha - beta   ;
    end
    
    % Résolution du problème quadratique
    LSk = [H jcdex' ; jcdex zeros(p, p)]        ;
    yk = -[ gfdex ; cdex] ;
    res=LSk\yk;
    dk = res(1:n);
    
    fdex_pred = fdex;
    gfdex_pred = gfdex;
    jcdex_pred = jcdex;
    tempx_pred = tempx;
    templambda_pred = templambda ;
    
    templambda = res(n+1:n+p);
    tempx = tempx +dk ;
    
    
    fdex= feval(une_f,tempx) ;
    cdex = feval(des_c,tempx) ;    
   
    for i =1:n
        gfdex(i, 1) =(feval(une_f, tempx +une_tol_h*E(:,i)) - fdex)/ une_tol_h;    
        jcdex(:, i) =(feval(des_c, tempx +une_tol_h*E(:,i)) - cdex) /une_tol_h;
    end
    
    % Update des valeurs sk_pred et yk_pred
    sk_pred = dk ;
    yk_pred = (gfdex + templambda * jcdex' ) - (gfdex_pred + templambda_pred * jcdex_pred');
    
    
    if (abs(dk) <une_tol_x)
        fin = 1;
    end

    if (norm(fdex - fdex_pred) < une_tol_f)
            fin = 2;
    end

    if (abs(gfdex) < une_tol_g)
        fin = 3;
    end
    if(k==un_nit_max)
        fin =  4                                                           ;
    end

    k = k +1 ;
end
x_opt    =   tempx                                                         ;
f_opt    =   fdex                                                      ;
nit      =   k                                                             ;
