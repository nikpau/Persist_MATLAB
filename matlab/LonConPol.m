classdef LonConPol
    
    properties
        
        % constant
        id
        maxPower = 1E6;
        upperSpeedBound = 8;
        maxDistToConsider = 300;
        vStrScale = 4;
        Amean = 1500;
        Ascale = 800;
        hMean = 5.5;
        hScale = 3;
        
        % observations
        obs
        obs_vStream
        obs_A
        obs_h
        distShipAhead
        wid;
    end
    
    methods
        function o = LonConPol(ID,ships,river)
            o.id = ID;
            o.obs_vStream = river.getMeanVStreamForVessel(ships, ID);
            o.wid = 0;
        end
        
        function o = computeObs(o,ships, river)
            ID = o.id;
            [shipAheadID, o.distShipAhead] = ships.findNextShipAhead(ID,1,ships.overtakeLevel(ID));
            
            
            if shipAheadID ~= 0
                vAhead = (ships.vx(shipAheadID));
                o.distShipAhead = o.distShipAhead - ships.length(ID)/2 - ships.length(shipAheadID)/2;

            else
                vAhead = o.upperSpeedBound;
                o.distShipAhead = o.maxDistToConsider;
            end
            o.distShipAhead = clip(o.distShipAhead,0,o.maxDistToConsider);
            
            % get river properties
            o.obs_vStream = river.getMeanVStreamForVessel(ships, ID);
            o.obs_A = river.getRiverProfileForVessel(ships,ID);
            o.obs_h = river.getWaterDepthForVessel(ships,ID);
            
            
            o.obs = [
                ships.dir(ID)* ships.vx(ID)/o.upperSpeedBound;
                o.distShipAhead /o.maxDistToConsider;
                ships.dir(ID)*(vAhead-ships.vx(ID))/o.upperSpeedBound;
                ships.P(ID)/o.maxPower;
                
                o.obs_vStream/o.vStrScale;
                (o.obs_A - o.Amean)/o.Ascale;
                (o.obs_h - o.hMean)/o.hScale;
                ];
        end
        
        function  [o,newAcc, newSquat, cf] = computeAccByPower(o,ships,ID,SimSet,river)
            
            obs_vStream_prev = o.obs_vStream;
            o.obs_vStream = river.getMeanVStreamForVessel(ships, ID);
            o.obs_A = river.getRiverProfileForVessel(ships, ID); 
            o.obs_h = river.getWaterDepthForVessel(ships,ID);
            
            % dynamische Variablen
            
                      % Breite Fluss/Kanal
            h=o.obs_h;            % Tiefe Fluss/Kanal
            A = o.obs_A;
            
            v = abs(ships.vx(ID));
            vRel = v - o.obs_vStream;
            
            % Konstanten
            rho=1000;       % Wasserdichte
            gamma=0.667;    % Voelligkeit/block coefficient
            Ls=ships.length(ID);         % GMS Laenge
            Bs=ships.width(ID);        % GMS Breite
            Ts=2.80;        % GMS "Abladetiefe" (Tiefgang)
            Tsmax=Ts;
            As = Bs*Ts;
            tSavg = 0.85; %%avg draught
            
            ms = 695258;
            
            ck = 0.8;
            ks=0.2;         % Rauigkeit Sohle [m]
            % (nicht angegeben, Wert so dass W3=W4 fuer vrel->0)
            
            kss = 0.0003;   % Unterschiff-Rauigkeit
            cStr = 0.78;
            alpha = 1.5E-4;
            
            
            
            n = A/As;
            vR = -v/(n - 1);
            mh = o.hydrodynMass(A,As,ms);
            
            sBa = 0.217*n^(-0.76)*gamma*vRel^2; % Squat nach Barrass
            newSquat = sBa;
            
            % Schiffswiderstand
            vReff = -vRel*(1/(A/(Bs*(Ts + sBa)) - 1));
            cW = 0.1+0.2*(Tsmax/h)^2;
            W1 = 0.5*cW*ck*rho*As*(vRel-vReff)^2;
            
            lambdaS = (1.89 + 1.62 * log10(h/kss))^(-2.5);
            W2 = 0.5*cW*ck*rho*lambdaS*Bs*Ls*(1+2*Ts/Bs)*gamma*(vRel-vReff)^2;
            
            Idir = -sign(o.obs_vStream-vReff);
            W3 = 0.023*Idir*cStr*rho*(ks*sBa/(4*h))^(1/3)*Ls*Bs*(Ts/h)*(1+2*Ts/Bs)*gamma*(o.obs_vStream-vReff)^2;
            
            W4 = sign(-vRel)*ms*9.81*alpha;
            
            W = W1 + W2 + W3 + W4;
            
            
            
            %
            % % power = 1000000;
            %      RHO = 1000;
            %
            % 	%%parameter water resistance
            % 	  c_a_cw = 1;
            % 	delta = 0.82;
            % 	c_K = 0.8;
            % 	c_Str = 0.78;
            % 	c_a_squat = 0.15;
            % 	c_vCr = 1;
            %
            % 	%%general parameter
            % 	bS = 9.5; %%width
            % 	mS = 695258; %%mass
            % 	lS = 105; %%vessel length
            % 	bSeff = 9.5; %%effective vessel width
            % 	tSmax = 1.4; %%max draught
            % 	tSavg = 0.85; %%avg draught
            %
            % 	mH = 75851; %%hydro mass
            
            %%thrust parameter
            thrustFactorA = 16.26;
            thrustFactorB = 6.25;
            thrustFactorC = 0.8;
            
       
            
            % river parameters
            %          aEP = 734.42; %% mittlere Querschnittsfläche des Ersatzprofils
            %          hEP = 4.23; %% avg water depth Ersatzprofil
            %          aSeff = 8.075; %% effektive frontale Querschnittsfläche
            %          bEP = 173.58; %% mittlere Breite des Ersatzprofils
            %          ks = 0.06; %%avg bottom friction
            %          kss = 3E-4; %%Hull roughness
            %          cG = 1;
            %          vvStr = 0; %% stream velocity over ground (schiffseigenes KS)
            %
            %
            %
            %         % compute water resistance
            %          cw = c_a_cw * (0.1 + 0.2 * (tSmax / hEP)^2);
            %          lambdaS = 4.0 * (1.89 - 1.62 * log10(kss / hEP))^(-2.5);
            %          vscr = 0.35 * c_vCr * sqrt(9.81 * hEP) * ((aEP / aSeff - 1.0))^0.25;
            %          ssq = (state.ships.vel(1) / vscr)^3 * c_a_squat * sqrt(tSavg * hEP);
            %
            %
            %          aEPa = aEP - ssq * bEP;
            %          vR = state.ships.vel(1) * (aEP / (aEPa - aSeff) - 1);
            %          vStrEff = vvStr - 1 * vR;
            %
            %
            %
            %          term1 = 0.5 * RHO * (state.ships.vel(1) + vR)^2 * bSeff * lS * (cw * tSmax / lS + 0.25 * lambdaS * (1.0 + 2.0 * tSavg / bSeff) * delta) * c_K;
            %          term2 = cG * 0.184 * (ks / (4.0 * hEP)^(1/3)) * RHO * 0.125 * vStrEff^2 * tSavg / hEP * lS * (bSeff + 2 * tSavg) * c_Str * delta;
            %          resistance = term1 + term2;
            
            
            
            %%compute thrust
            thrust = thrustFactorC * (thrustFactorA * ships.P(ID)) / (ships.P(ID)^(1/3) + thrustFactorB *vRel);
            

            A_next = 1000;
            A_prev = 1000;
            %%compute acceleration
            dvStr_dt = (o.obs_vStream - obs_vStream_prev)/SimSet.dT;
            
            mhPrev = o.hydrodynMass(A_prev,As,ms);
            mhNext = o.hydrodynMass( A_next,As,ms);
            dmh_dt = (mhNext-mhPrev)/SimSet.dT;
            
            % Bewegungsgleichung
            newAcc = 1.0 / (ms + mh) * (thrust - W - vRel * dmh_dt + mh * dvStr_dt);
            o.wid = newAcc;
            
                 %%Cf parameter
            	gamma = 0.82;
            	gammaL = 0.8;
            	c_ast = 3.2;
            	c_i = 0.5;
            	f_v2 = 0.3;
            	c_korr = 1.4;
            	beta = 2.0;
            	betaC = 0.5;
            	delta_c_fs = 0;
            	delta_c_f_korr_BS = 0;
            
            
            % compute CF value
   
         cD0 = 0.22 * (Ls / (gamma * Bs))^betaC;
         cDmax = cD0 + (c_ast - cD0) * (tSavg / h)^beta;
         vStr_vSdW = o.obs_vStream  / vRel;
         f1 = (1.0 + 2.0 * vStr_vSdW) + f_v2 * (vStr_vSdW)^2;

         f11 = (tSavg / h) / (1.0 - (tSavg / h));
         f2 = 1.0 + 0.35 * log(Bs / tSavg + 1.0);
         f3 = 0.5 * pi * tSavg^2 / (tSavg * Bs + 0.2 * tSavg * tSavg);
         f4 = (tSavg * Bs + 0.2 * tSavg^2) / (tSavg * Bs);
         cmhy = f4 * (f3 * f2 + f11);

         kPrime = 1.0 - (3.4^(-0.22 * (Ls / Bs - 1.0)));
         cf = 0.5 + c_i * (Bs / Ls) + (2.0 / cDmax * Bs / Ls * gamma / gammaL * c_korr * f1 + 1.0 / 6.0) / (gamma / gammaL * Bs / Ls * cmhy / cDmax * 4 * kPrime + 1.0);
        cf = cf + delta_c_fs + delta_c_f_korr_BS;
            
        end
        
        function mh = hydrodynMass(o,A,As,ms)
            
            n = A/As;
            mh = ms*(1/(n-1) + 0.1);
            
        end
        
    end
end

