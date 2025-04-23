classdef pche < matlab.mixin.Copyable
    %RECALCULATION HEAT EXCHANGER OF TYPE PRINTED CIRCUIT HEAT EXCHANGER
    %references: VDI Waermeatlas, 11.Auflage, SpringerVieweg
                  
     % charge/discharge
     %1: process medium 1: index 1, hot side (CO2, water, or air)
     %2: Steel
     %3: process medium 2: index 3, cold side (CO2)
     
    properties
        calculated = 0
        name
        flow_type
        process_medium1
        process_medium3
        prop1;
        prop3;

        HTcorr;      % heat transfer correlation; available: "overallHTC", "meshram", "meshram2", "kim", "saeedkim", "cheng", "dittusboelter", "jacksonpitla", "gnielinski_filonenko", "gnielinski_konakov"
        % references
        % Meshram et al. http://dx.doi.org/10.1016/j.applthermaleng.2016.05.033
        % Kim et al. http://dx.doi.org/10.1016/j.anucene.2016.01.019
        % Saeed Kim https://doi.org/10.1016/j.enconman.2019.04.058
        % Cheng et al. https://doi.org/10.1016/j.applthermaleng.2021.116882
        % Jackson http://dx.doi.org/10.1016/j.nucengdes.2012.09.040; 
        % Pitla, Groll, Ramadhyani 2002 New correlation to predict the heat
        % transfer coefficient during in-tube cooling of turbulent supercritical CO2, Int. Journal of Refrigeration
        % Gnielinski VDI Heat Atlas 
        
        QSum;
        % start conditions
        Qdot_start = 0;
        T1_start = 0;
        T3_start = 0;

        mDot1;       % mass flow process medium 1(kg/s)               [1,1]
        mDot3;       % mass flow process medium 2(kg/s)               [1,1]

        d1;          % diameter of channel area (zigzag) (m)
        w1=0;        % distance between quarter circles for "other_zigzag"; for semicircular: w=0
        pitch1;      % zigzag channel pitch (m)
        alpha1;      % zigzag channel angle (m)   
        dh1;         % hydraulic diameter (m)        
        nChannel1;   % number of channels per plate
        nPlates1;    % number of plates
        ePlate1;     % total thickness of plate (m)
        d3;          % diameter of channel area (m)
        w3=0;        % distance between quarter circles for "other_zigzag"; for semicircular: w=0
        pitch3;      % zigzag channel pitch
        alpha3;      % zigzag channel angle   
        dh3;         % hydraulic diameter        
        nChannel3;   % number of channels per plate
        nPlates3;    % number of plates
        ePlate3;     % total thickness of plate
        l;           % effective hex length (m)                       [1,1]
        l_in;        % inlet/outlet length as a fraction of total heat exchanger length where cold channels change direction and make a cross-counter heat exchanger (per inlet and outlet)
        nCells;      % number of cells (-)                            [1,1]

        T1left;      % temperature left side (K)                      [1,n]
        T1right;     % temperature right side (K)                     [1,n]
        T1           % average temperature (K)                        [1,n]
        T3left;      % temperature left side (K)                      [1,n]
        T3right;     % temperature right side (K)                     [1,n]
        T3           % average temperature (K)                        [1,n]
        rho1;        % density (kg/m^3)                               [1,n]
        rho3;        % density (kg/m^3)                               [1,n]
        p1           % average pressure (Pa)                          [1,n]
        delta_p1;    % pressure loss (Pa)                             [1,n]
        p1in         % initial pressure (Pa)                          [1,1]
        h1left;      % specific enthalpy left side (J/kg)             [1,n]
        h1right;     % specific enthalpy right side (J/kg)            [1,n]
        p3           % average pressure (Pa)                          [1,n]
        delta_p3;    % pressure loss (Pa)                             [1,n]
        p3in;        % initial pressure (Pa)                          [1,1]
        h3left;      % specific enthalpy left side (J/kg)             [1,n]
        h3right;     % specific enthalpy right side (J/kg)            [1,n]
        k            % overall heat transfer coefficient(W/m^2*K)     [1,n]        
        al1;         % heat transfer coefficient h  (W/m^2*K)         [1,n]
        al3;         % heat transfer coefficient h  (W/m^2*K)         [1,n]
        Re1;         % Re number                                      [1,n]
        Re3;         % Re number                                      [1,n]
        Pr1;         % Pr number                                      [1,n]
        Pr3;         % Pr number                                      [1,n]

        Qdot1;       % heat flow (W)                                [1,n+1]
        Qdot_1       % heat flow medium 1 (W)                         [1,1]
        Qdot_3       % heat flow medium 2(W)                          [1,1]
        QLoss_Rec;   % heat flow loss (W)                             [1,n]
        F=1;         % flow characteristics crosscurrent/countercurrent [1,1]                                                                                                          
        lambda_st;   % heat conductivity for steel tube(W/m*K)        [1,1]
        charge=0;    % logical value (-)                              [1,n]
        kA;
        Atot;

        calctime;      % elapsed time
        count_it;      % count of iterations
        checkRe1;      % range of validity for dimensionless numbers met?
        checkRe3;
        checkPr1;
        checkPr3; 
        checkT1;
        checkT3; 
    end
    
    methods
        function obj = pche(prop1,prop3,mDot1,mDot3)
            if nargin==0
                return
            else
                obj.prop1 = prop1;
                obj.prop3 = prop3;
                obj.mDot1 = mDot1;
                obj.mDot3 = mDot3;
            end

        end
        
        %% TQ DIAGRAM
        % plot TQ-diagram
        function TQ_diagram(obj)
            figure(1)
            if obj.process_medium1 == "CO2"
                obj.prop1 = CO2();
            elseif obj.process_medium1 == "water" % not in use. following code not equipped for water
                obj.prop1 = IF97();
            end

            if obj.process_medium3 == "CO2"
                obj.prop3 = CO2();
            elseif obj.process_medium3 == "water" % not in use. following code not equipped for water
                obj.prop3 = IF97();
            end

            len = sum([obj.nCells]);
            % a = T1in ~= 0;
            % T1in = T1in(a);
            % p = [obj.p1in];
            % p = p(a);
            rho1in = obj.prop1.rho(obj.p1in, obj.T1_start);
            % b = T3in ~= 0;
            % T3in = T3in(b);            

            if obj.T1_start > obj.T3_start  %charge 
                [obj.charge]=deal(true); 
                Qdot_start = [0, obj.Qdot_start];
            else           %discharge
                [obj.charge]=deal(false);
                Qdot_start = [0, -obj.Qdot_start];
   
            end

            % pressures
            if size(obj.delta_p1,1) == 1 && size(obj.delta_p1,2) == 1
                obj.p1=linspace(obj.p1in,obj.p1in-obj.delta_p1,obj.nCells);
            else
                obj.p1=[obj.p1in,obj.p1in*ones(1,len)-cumsum(obj.delta_p1)];
            end

            if size(obj.delta_p3,1) == 1 && size(obj.delta_p3,2) == 1
                obj.p3=linspace(obj.p3in,obj.p3in-obj.delta_p3,obj.nCells);
            else
                obj.p3=[obj.p3in,obj.p3in*ones(1,len)-cumsum(obj.delta_p3)];
            end
            
            % process medium 2
            QDot3=linspace(0,Qdot_start(end),obj.nCells+1);

            H3In=obj.prop3.h(obj.p3in,obj.T3_start); % enthalpy
            if obj.charge
                obj.h3right=H3In+fliplr(QDot3(1:end-1))./obj.mDot3;
                obj.h3left=H3In+fliplr(QDot3(2:end))./obj.mDot3;
            else
                obj.h3left=H3In+(QDot3(1:end-1))./obj.mDot3;
                obj.h3right=H3In+(QDot3(2:end))./obj.mDot3;
            end
            
            %
            obj.p3;
            obj.T3left=obj.prop3.T_ph(obj.p3,obj.h3left); % temperature
            obj.T3right=obj.prop3.T_ph(obj.p3,obj.h3right);
            obj.T3=(obj.T3left + obj.T3right)./2;

            % process medium 1  
            obj.h1left = NaN(1,len);
            obj.h1right = NaN(1,len);
            H1In=obj.prop1.h_rhoT(rho1in,obj.T1_start); % enthalpy

            QDot1m=linspace(0,Qdot_start(end),obj.nCells+1);
            if ~obj.charge
                QDot1m = fliplr(QDot1m);
            end
              
            obj.h1left=H1In-QDot1m(1:end-1)./obj.mDot1;
            obj.h1right=H1In-QDot1m(2:end)./obj.mDot1;
            obj.T1left=obj.prop1.T_ph(obj.p1,obj.h1left); % temperature
            obj.T1right=obj.prop1.T_ph(obj.p1,obj.h1right);
               

            obj.T1=(obj.T1left+obj.T1right)./2;

            if obj.process_medium1 == "CO2" && obj.process_medium3 == "CO2"
                obj.rho1 = obj.prop1.rho(obj.p1,obj.T1); % for CO2
                obj.rho3 = obj.prop3.rho(obj.p3,obj.T3); % for CO2
            end

           % plot 
            if obj.T1_start > obj.T3_start % charge
                p_1 = plot(abs(QDot3),[[obj.T1left-273.15,obj.T1right(end)-273.15];[obj.T3left(1)-273.15,obj.T3right-273.15]]);
            else            % discharge
                p_1 = plot(abs(QDot3),[([obj.T1left(1)-273.15,obj.T1right-273.15]);([obj.T3left-273.15,obj.T3right(end)-273.15])]);
            end
            
            if obj.process_medium1 == "CO2" && obj.process_medium3 == "CO2"
                legend(p_1,{'CO_2', 'CO_2'})
            end

            xlabel("QDot in W")
            ylabel("T in °C")

            obj.Qdot1 = QDot3;
                         
        end

        %% HEAT EXCHANGER CALCULATION
        % DESCRIPTIVE TEXT 
        function hex_calculation(obj,idx)
            spec1 = obj(1).prop1;
            spec3 = obj(1).prop3;

            % get values
            len = sum([obj.nCells]);
            Rho1 = [obj.rho1];
            Rho3 = [obj.rho3];
            P1 = [obj.p1];
            P3 = [obj.p3];
            T1m = [obj.T1];
            T3m = [obj.T3];
            T1in=[obj.T1left];
            T3in=[obj.T3right];
            h1in=[obj.h1right];
            h3in=[obj.h3left];

            if obj(1).flow_type == "crosscurrent"
                obj(1).F = 0.869;
            elseif obj(1).flow_type == "countercurrent"
                obj(1).F = 1;
            elseif obj(1).flow_type == "custom"
                obj(1).F = ones(1,len);
                inl=round(obj.l_in*obj.nCells);
                obj(1).F(1:inl)=0.869-(1-0.869)/(inl-1)+(1-0.869)/(inl-1).*(1:inl); 
                obj(1).F(end-inl+1:end)=flip(obj(1).F(1:inl));                      
            end

            mdot1 = obj.mDot1; 
            mdot3 = obj.mDot3;
            HT = obj.HTcorr;
            s = (obj.ePlate1+obj.ePlate3)/2;

            % geometric data & calculations
            Area1=obj.d1^2*pi()/8+obj.d1/2*obj.w1;
            Peri1=obj.d1*pi()/2+obj.d1+2*obj.w1;
            obj.dh1=4*Area1/Peri1;                  % hydraulic diameter

            Area3=obj.d3^2*pi()/8+obj.d3/2*obj.w3;
            Peri3=obj.d3*pi()/2+obj.d3+2*obj.w3;
            obj.dh3=4*Area3/Peri3;                  % hydraulic diameter

            A1=Peri1.*obj.nChannel1*obj.nPlates1.*obj.l./obj.nCells;  % hex area over all channels per cell, side 1
            A3=Peri3.*obj.nChannel3*obj.nPlates3.*obj.l./obj.nCells;  % hex area over all channels per cell, side 3
            Am=(A1-A3)./log(A1./A3);      % area for heat conduction calculation
  
            if ~obj(1).charge %discharge
                T1in=fliplr(T1in);
                T3in=fliplr(T3in);
                h1in=fliplr(h1in);
                h3in=fliplr(h3in);      
                T1m=fliplr(T1m);
                T3m=fliplr(T3m);
                P1=fliplr(P1);
                Rho1=fliplr(Rho1);
                mdot1=fliplr(mdot1);
            end
         
            err=repmat(10,1,5);
            QDot=0;
            count=0;
            isConv=false;
            disp([obj.name]);
            QDoti=-(obj.Qdot1(1:end-1)-obj.Qdot1(2:end));
            al1=3000;
            al3=3000;
           
            while err(end)>1 || ~isConv
                
                if HT == "meshram2" || HT == "meshram" || HT == "kim" || HT=="saeedkim" || HT=="cheng" || HT == "overallHTC" || HT == "gnielinski_filonenko" || HT == "gnielinksi_konakov" || HT == "dittusboelter"
                    %dimensionless numbers
                    Rho1=spec1.rho(P1,T1m);
                    Rho3=spec3.rho(P3,T3m);

                    eta1 = spec1.my(Rho1,T1m);          % dynamic viscosity Pa*s
                    lambda1 = spec1.lambda(Rho1,T1m);   % heat conductivity
                    u1 = mdot1./(Area1.*obj.nChannel1*obj.nPlates1)./Rho1;      % velocity
                    Pr1 = spec1.Pr(Rho1,T1m);           % Pr-number
                    Re1 = u1.*obj.dh1.*Rho1./eta1;      % Re-number
        
                    eta3 = spec1.my(Rho3,T3m);          % viscosity
                    lambda3 = spec1.lambda(Rho3,T3m);   % heat conductivity
                    u3 = mdot3./(Area3.*obj.nChannel3*obj.nPlates3)./Rho3;      % velocity
                    Pr3 = spec3.Pr(Rho3,T3m);           % Pr-number
                    Re3 = u3.*obj.dh3.*Rho3./eta3;      % Re-number

                elseif HT == "jacksonpitla"
                    Rho1=spec1.rho(P1,T1m);
                    Rho3=spec3.rho(P3,T3m);

                    T3mw=T3m+QDoti./(al3.*A3.*obj.F);                    
                    T1mw=T1m-QDoti./(al1.*A1.*obj.F);
                    meanT=mean([T1mw; T3mw],1);
                    T3mw(T1mw<T3mw)=meanT(T1mw<T3mw);
                    T1mw(T1mw<T3mw)=meanT(T1mw<T3mw);

                    %Jackson
                    Tpc3=CO2.Tpc(P3);

                    idx3=NaN(size(T3m));
                    n_jack=NaN(size(T3m));

                    idx3(T3m<T3mw & T3mw<Tpc3)=1;
                    idx3(1.2*Tpc3<T3m & T3m<T3mw)=1;
                    idx3(T3m<Tpc3 & Tpc3<T3mw)=2;
                    idx3(Tpc3<T3m & T3m<1.2*Tpc3 & T3m<T3mw)=3;
                    
                    n_jack(idx3==1)=0.4;
                    n_jack(idx3==2)=0.4+0.2.*(T3mw(idx3==2)./Tpc3(idx3==2)-1);
                    n_jack(idx3==3)=0.4+0.2.*(T3mw(idx3==3)./Tpc3(idx3==3)-1).*(1-5.*(T3m(idx3==3)./Tpc3(idx3==3)-1));
                    
                    Cpb= spec3.cp(P3, T3m);
                    hb=spec3.h(P3,T3m);
                    hw=spec3.h(P3,T3mw);
                    Cpm=(hw-hb)/(T3mw -T3m);
                    
                    u3 = mdot3./(Area3.*obj.nChannel3*obj.nPlates3)./Rho3;  % velocity
                    eta3 = spec1.my(Rho3,T3m);      % viscosity
                    Pr3 = spec3.Pr(Rho3,T3m);       % Pr-number
                    Re3 = u3.*obj.dh3.*Rho3./eta3;  % Re-number
                    
                    % Pitla/Gnielinski
                    T1f=(T1m+T1mw)./2;              % film temperature

                    Rho1f = spec1.rho(P1,T1f);
                    Rho1w = spec1.rho(P1,T1mw);

                    eta1 = spec1.my(Rho1f,T1m);     % dynamic viscosity Pa*s
                    u1 = mdot1./(Area1.*obj.nChannel1*obj.nPlates1)./Rho1;  % velocity
                    Pr1 = spec1.Pr(Rho1,T1m);       % Pr-number
                    Re1 = u1.*obj.dh1.*Rho1./eta1;  % Re-number

                    eta1w = spec1.my(Rho1w,T1mw);   % dynamic viscosity Pa*s
                    u1w = mdot1./(Area1.*obj.nChannel1*obj.nPlates1)./Rho1w; % velocity
                    Pr1w = spec1.Pr(Rho1,T1mw);     % Pr-number
                    Re1w = u1w.*obj.dh1.*Rho1w./eta1w; % Re-number

                    lambda1 = spec1.lambda(Rho1,T1m); % heat conductivity
                    lambda3 = spec1.lambda(Rho3,T3m); % heat conductivity  

                elseif HT == "dittusboelter_film"
                    Rho1=spec1.rho(P1,T1m);
                    Rho3=spec3.rho(P3,T3m);

                    T1mw=T1m-QDoti./(al1.*A1.*obj.F);
                    T3mw=T3m+QDoti./(al3.*A3.*obj.F);
                    T1f=(T1m+T1mw)./2;                 % film temperature
                    T3f=(T3m+T3mw)./2;
                    
                    Rho1f = spec1.rho(P1,T1f);
                    Rho3f = spec1.rho(P3,T3f);

                    eta1 = spec1.my(Rho1f,T1f);         % dynamic viscosity Pa*s
                    lambda1 = spec1.lambda(Rho1f,T1f);  % heat conductivity
                    u1 = mdot1./(Area1.*obj.nChannel1*obj.nPlates1)./Rho1f;     % velocity
                    Pr1 = spec1.Pr(Rho1f,T1f);          % Pr-number
                    Re1 = u1.*obj.dh1.*Rho1f./eta1;     % Re-number
    
                    eta3 = spec1.my(Rho3f,T3f);         % viscosity
                    lambda3 = spec1.lambda(Rho3f,T3f);  % heat conductivity
                    u3 = mdot3./(Area3.*obj.nChannel3*obj.nPlates3)./Rho3f;     % velocity
                    Pr3 = spec3.Pr(Rho3f,T3f);          % Pr-number
                    Re3 = u3.*obj.dh3.*Rho3./eta3;      % Re-number
                    
                else
                    error('No HT corr fitting to input');
                end

                if HT == "meshram2"
                    [al1,al3,f1,f3]=meshram2();

                elseif HT == "meshram"
                    [al1,al3]=meshram();

                elseif HT == "kim"
                    [al1,al3]=kim();

                elseif HT == "saeedkim"
                    [al1,al3]=saeedkim();

                elseif HT == "cheng"
                    [al1,al3]=cheng();

                elseif HT == "dittusboelter" || HT == "dittusboelter_film"
                    [al1,al3]=dittusboelter();

                elseif HT == "gnielinski_filonenko"
                    %friction factor according to Filonenko
                    f1=1./(1.82.*log10(Re1)-1.64).^2;
                    f3=1./(1.82.*log10(Re3)-1.64).^2;
                    %Darcy friction factor
                    % f1=1./(0.79.*log10(Re1)-1.64).^2;
                    % f3=1./(0.79.*log10(Re3)-1.64).^2;
                    [al1,al3]=gnielinski();

                elseif HT == "gnielinksi_konakov"
                    %friction factor according to Konakov
                    f1=1./(1.8.*log10(Re1)-1.5).^2;
                    f3=1./(1.8.*log10(Re3)-1.5).^2;                    
                    [al1,al3]=gnielinski();

                elseif HT == "jacksonpitla"
                    [al1,al3]=jacksonpitla();
                
                elseif HT == "overallHTC"
                    [al1,al3]=overallHTC();

                end

                dT = (T1m - T3m);

                k=((1./(al1.*A1)) + (s./(obj.lambda_st.*Am))+ 1./(al3.*A3)).^-1;

                % transferred heat per cell
                QDoti=obj.F.*k.*dT;
                               
                QDot3=QDoti;

                %energy balance 
                h1out=h1in-QDoti./mdot1;
                h1in(2:end)=h1out(1:end-1);
                
                T1out=spec1.T_ph(P1,h1out); 
                               
                %Consistency check
                if obj(1).charge
                    check=T1out<T3in;
                else
                    check=T1out>T3in;    
                end
                T1out(check)=T3in(check);
                
                T1in(2:end)=T1out(1:end-1);
                T1m=mean([T1in;T1out],1);
               
                h3out=h3in+QDot3./mdot3;
                h3in(1:end-1)=h3out(2:end);
                
                T3out=spec3.T_ph(P3,h3out);
                
                %Consistency check
                if obj(1).charge
                    check=T3out>T1in;
                else
                    check=T3out<T1in;    
                end
                T3out(check)=T3out(check);
                
                T3in(1:end-1)=T3out(2:end);
                T3m=mean([T3in;T3out],1);
               
                % pressure loss
                delta_p1 = 2.*Rho1.*obj.l/obj.nCells.*u1.^2./obj.dh1; % FIX! deltap_i();
                P1 = obj.p1in*ones(1,len); %FIX P1 = obj.p1in*ones(1,len) -  delta_p1;
                delta_p3=0;
                P3 = obj.p3in*ones(1,len) -  delta_p3;
               

                if obj(1).process_medium1 == "water"
                    Rho1 = hex.rho(P1,mean([h1in;h1out])); %for H2O
                elseif obj(1).process_medium1 == "CO2"
                    Rho1 = spec1.rho(P1,T1m); % for CO2
                end
                Rho3 = spec3.rho(P3,T3m);
              
                % termination condition
                QDotsum=sum(QDot3);
              
                err=circshift(err,-1);
                err(end)=abs(QDotsum-QDot);
                QDot=QDotsum;
           
                % output information
                if mod(count,10)==0
                   disp([obj.name+": count="+count+", err="+err(end)+" QDot="+QDot])
                end
                
                % termination if 10000 iterations
                count=count+1;
                if count>10000
                    break
                end

                %convergence condition
                derr=diff(err);
                isConv=(all(derr<0) && err(end)/err(1)>0.7) || max(abs(derr))<0.2;
            end 
            % end of iteration
            %% save values
            
            if ~obj(1).charge
                T1m = fliplr(T1m);
                T3m = fliplr(T3m);
                T1in = fliplr(T1in);
                T1out = fliplr(T1out);
                T3out = fliplr(T3out);
                T3in = fliplr(T3in);
                h1in = fliplr(h1in);
                h1out = fliplr(h1out);
                h3out = fliplr(h3out);
                h3in = fliplr(h3in);
                P1 = fliplr(P1);
                delta_p1 = fliplr(delta_p1);
                QDoti = fliplr(QDoti);
                QDot3 = fliplr(QDot3);
            end
            step = 0;
            
            for t = 1:size(obj,2)
                if obj(1).charge
                    obj(t).T1 = T1m(step+1:step+obj(t).nCells);
                    obj(t).T3 = T3m(step+1:step+obj(t).nCells);

                    obj(t).T1left = T1in(step+1:step+obj(t).nCells);
                    obj(t).T1right = T1out(step+1:step+obj(t).nCells);
                    obj(t).T3left = T3out(step+1:step+obj(t).nCells);
                    obj(t).T3right = T3in(step+1:step+obj(t).nCells);

                    obj(t).h1left = h1in(step+1:step+obj(t).nCells);
                    obj(t).h1right = h1out(step+1:step+obj(t).nCells);
                    obj(t).h3left = h3out(step+1:step+obj(t).nCells);
                    obj(t).h3right = h3in(step+1:step+obj(t).nCells);
                else
                    obj(t).T1 = T1m(step+1:step+obj(t).nCells);
                    obj(t).T3 = T3m(step+1:step+obj(t).nCells);

                    obj(t).T1left = T1out(step+1:step+obj(t).nCells);
                    obj(t).T1right = T1in(step+1:step+obj(t).nCells);
                    obj(t).T3left = T3in(step+1:step+obj(t).nCells);
                    obj(t).T3right = T3out(step+1:step+obj(t).nCells);

                    obj(t).h1left = h1out(step+1:step+obj(t).nCells);
                    obj(t).h1right = h1in(step+1:step+obj(t).nCells);
                    obj(t).h3left = h3in(step+1:step+obj(t).nCells);
                    obj(t).h3right = h3out(step+1:step+obj(t).nCells);
                end

                obj(t).p1 = P1(step+1:step+obj(t).nCells);
                obj(t).delta_p1 = delta_p1;                    
                obj(t).Qdot_1 = QDoti(step+1:step+obj(t).nCells);
                obj(t).Qdot_3 = QDot3(step+1:step+obj(t).nCells);
                obj(t).k = k;
                obj(t).rho1=Rho1;
                obj(t).rho3=Rho3;
                obj(t).al1=al1;
                obj(t).al3=al3;
                obj(t).Re1=Re1;
                obj(t).Re3=Re3;
                obj(t).QSum = sum(QDoti(step+1:step+obj(t).nCells));
                obj(t).count_it=count;
                obj(t).calculated=1;

                step = step + obj(t).nCells;
            end
             
           
            %% HEAT TRANSFER CORRELATIONS
            function [alpha_hot, alpha_cold]=overallHTC()
                alpha_hot = ones(1,len);    %index 1
                alpha_cold = ones(1,len);   %index 3
                alpha_hot = 1398.*alpha_hot;
                alpha_cold = 1260.5.*alpha_cold;
            end

            function [alpha_hot, alpha_cold]=meshram()
                Nu_hot = 87.56*(obj.pitch1/0.012)^-0.178*(obj.alpha1/116)^-0.9306*ones(1,len); 
                Nu_cold = 85.95*(obj.pitch3/0.012)^-0.171*(obj.alpha3/116)^-0.8912*ones(1,len); 

                alpha_hot = Nu_hot .* lambda1 ./ obj.dh1;
                alpha_cold = Nu_cold .* lambda3 ./ obj.dh3;

                obj.checkT1=T1m<630 & T1m>470;
                obj.checkRe1=Re1<32000 & Re1>5000;
                obj.checkT3=T3m<520 & T3m>400;
                obj.checkRe3=Re3<32000 & Re3>5000;
            end

            function [alpha_hot, alpha_cold, frictionfactor_hot, frictionfactor_cold]=meshram2()
                % Nu = a * Re^b * Pr^c
                Nu1_hot = 0.0174 .* Re1.^0.893 .* Pr1.^0.7; %valid at 470<T1m<630 
                Nu2_hot = 0.0205 .* Re1.^0.869 .* Pr1.^0.7; %valid at 580<T1m<730 
                Nu1_cold = 0.0177 .* Re3.^0.871 .* Pr3.^0.7; %valid at 400<T3m<520 
                Nu2_cold = 0.0213 .* Re3.^0.876 .* Pr3.^0.7; %valid at 500<T3m<640 

                isT1low=T1m<580;
                isT1middle=T1m>580&T1m<630;
                isT1high=T1m>630;
                Nu_hot(isT1low)=Nu1_hot(isT1low);
                Nu_hot(isT1middle)=(Nu1_hot(isT1middle)+Nu2_hot(isT1middle))./2;
                Nu_hot(isT1high)=Nu2_hot(isT1high);

                isT3low=T3m<500;
                isT3middle=T3m>500&T3m<520;
                isT3high=T3m>520;
                Nu_cold(isT3low)=Nu1_cold(isT3low);
                Nu_cold(isT3middle)=(Nu1_cold(isT3middle)+Nu2_cold(isT3middle))./2;
                Nu_cold(isT3high)=Nu2_cold(isT3high);
                
                alpha_hot = Nu_hot .* lambda1 ./ obj.dh1;
                alpha_cold = Nu_cold .* lambda3 ./ obj.dh3;
                
                % f = f(Re)
                frictionfactor1_hot=0.867.*Re1.^(-0.522)+0.04;
                frictionfactor2_hot=0.819.*Re1.^(-0.671)+0.044;
                frictionfactor1_cold=0.869.*Re1.^(-0.512)+0.041;
                frictionfactor2_cold=0.804.*Re1.^(-0.711)+0.045;                
                
                frictionfactor_hot(isT1low)=frictionfactor1_hot(isT1low);
                frictionfactor_hot(isT1middle)=(frictionfactor1_hot(isT1middle)+frictionfactor2_hot(isT1middle))./2;
                frictionfactor_hot(isT1high)=frictionfactor2_hot(isT1high);        

                frictionfactor_cold(isT1low)=frictionfactor1_cold(isT1low);
                frictionfactor_cold(isT1middle)=(frictionfactor1_cold(isT1middle)+frictionfactor2_cold(isT1middle))./2;
                frictionfactor_cold(isT1high)=frictionfactor2_cold(isT1high); 
                
                obj.checkT1=T1m<630 & T1m>470;
                obj.checkRe1=Re1<32000 & Re1>5000;
                obj.checkT3=T3m<520 & T3m>400;
                obj.checkRe3=Re3<32000 & Re3>5000;
            end

            function [alpha_hot, alpha_cold]= kim()
                % Nu = a * Re^b
                Nu_hot = 0.0292 .* Re1.^0.8138; 
                Nu_cold = 0.0188 .* Re3.^0.8742; 

                alpha_hot = Nu_hot .* lambda1 ./ obj.dh1;
                alpha_cold = Nu_cold .* lambda3 ./ obj.dh3;
                
                obj.checkRe1=Re1<58000 & Re1>2000;
                obj.checkRe3=Re3<55000 & Re3>2000;
                obj.checkPr1=Pr1<1 & Pr1>0.7;
                obj.checkPr3=Pr3<1 & Pr3>0.7;             
            end

            function [alpha_hot, alpha_cold]= saeedkim()
                % Nu = a * Re^b
                Nu_hot = 0.041 .* Re1.^0.83.* Pr1.^0.95; 
                Nu_cold = 0.041 .* Re3.^0.83.* Pr3.^0.95;

                alpha_hot = Nu_hot .* lambda1 ./ obj.dh1;
                alpha_cold = Nu_cold .* lambda3 ./ obj.dh3;
                
                obj.checkRe1=Re1<60000 & Re1>3000;
                obj.checkRe3=Re3<60000 & Re3>3000;
                obj.checkPr1=Pr1<1.2 & Pr1>0.7;
                obj.checkPr3=Pr3<1.2 & Pr3>0.7;             
            end

            function [alpha_hot, alpha_cold]= cheng()
                % Nu = a * Re^b
                Nu_hot = 0.02475 .* Re1.^0.7621; 
                Nu_cold = 0.02063 .* Re3.^0.7678;

                alpha_hot = Nu_hot .* lambda1 ./ obj.dh1;
                alpha_cold = Nu_cold .* lambda3 ./ obj.dh3;
                
                obj.checkRe1=Re1<23888 & Re1>4897;
                obj.checkRe3=Re3<15631 & Re3>3213;
                obj.checkPr1=Pr1<0.784 & Pr1>0.765;
                obj.checkPr3=Pr3<1.1 & Pr3>1.01;             
            end

            function [alpha_hot, alpha_cold]=gnielinski() 
                Nu_hot=(f1./8.*(Re1-1000).*Pr1./(1+12.7.*(Pr1.^(2/3)-1).*sqrt(f1./8))); 
                Nu_cold=(f3./8.*(Re3-1000).*Pr3./(1+12.7.*(Pr3.^(2/3)-1).*sqrt(f3./8))); 

                alpha_hot = Nu_hot .* lambda1 ./ obj.dh1;
                alpha_cold = Nu_cold .* lambda3 ./ obj.dh3;

                obj.checkRe1=Re1<5e6 & Re1>3000;
                obj.checkRe3=Re3<5e6 & Re3>3000;
                obj.checkPr1=Pr1<2000 & Pr1>0.5;
                obj.checkPr3=Pr3<2000 & Pr3>0.5; 

                %f_c=(1/(1.8*log10(Re)-1.5)).^2;            
            end

            function [alpha_hot, alpha_cold]= dittusboelter()
                % Nu = a * Re^b
                Nu_hot = 0.023 .* Re1.^0.8.*Pr1.^0.4;
                Nu_cold = 0.023 .* Re3.^0.8.*Pr3.^0.4;

                alpha_hot = Nu_hot .* lambda1 ./ obj.dh1;
                alpha_cold = Nu_cold .* lambda3 ./ obj.dh3;
            end

            function [alpha_hot, alpha_cold]= jacksonpitla()
                % Jackson
                rhob=spec3.rho(P3,T3m);
                rhow=spec3.rho(P3,T3mw);
                Nu_cold=0.0183.*Re3.^0.82.*Pr3.^0.5+(rhow./rhob).^0.3.*(Cpm./Cpb).^n_jack;
                
                % Pitla (using Gnielinksi correlation and Filonenko friction factor)
                f1=1./(0.79.*log10(Re1)-1.64).^2;
                Nu_b=(f1./8.*(Re1-1000).*Pr1./(1+12.7.*(Pr1.^(2/3)-1).*sqrt(f1./8)));
                f1w=1./(0.79.*log10(Re1w)-1.64).^2;
                Nu_w=(f1w./8.*(Re1w-1000).*Pr1w./(1+12.7.*(Pr1w.^(2/3)-1).*sqrt(f1w./8)));
                kw=CO2.lambda(Rho1w,T1mw);
                kb=CO2.lambda(Rho1,T1m);
                Nu_hot=((Nu_w+Nu_b)./2).*kw./kb;

                alpha_hot = Nu_hot .* lambda1 ./ obj.dh1;
                alpha_cold = Nu_cold .* lambda3 ./ obj.dh3;
            end


        end

        function plot_Tx(obj,type)
            if nargin < 2
                type = 1;
            end
            T1l = [obj.T1left];
            T1r = [obj.T1right];
            T3l = [obj.T3left];
            T3r = [obj.T3right];
            step = 0;
            step2 = 0;
            XI = NaN(1,sum([obj.nCells]));
            for n = 1:size(obj,2)
                XI(step+1:step+obj(n).nCells) = step2 + linspace(obj(n).l./obj(n).nCells,obj(n).l,obj(n).nCells);
                if type == 2
                    p_1 = plot([step2 XI(step+1:step+obj(n).nCells)], [T1l(step+1) T1r(step+1:step+obj(n).nCells)]-273.15,'Color',[0 0.4470 0.7410]);
                    hold on
                end
                step = step + obj(n).nCells;
                step2 = XI(step);
            end
            XII = [0 XI];
            if type == 1
                p_1 = plot(XII,[T1l(1) T1r]-273.15,'Color',[0 0.4470 0.7410]);
                hold on
            end
            p_3 = plot(XII,[T3l(1) T3r]-273.15,'Color',[0.9290 0.6940 0.1250]);
            hold off
            xlabel("x in m")
            ylabel("T in C")
        end

        function plot_TQ(obj, type)
            if nargin < 1
                type = 1;
            end
            Q = cumsum(abs([obj.Qdot_1]*10^-6));
            T1l = [obj.T1left];
            T1r = [obj.T1right];
            T3l = [obj.T3left];
            T3r = [obj.T3right];   
            QII = [0 Q];

            step = 0;
            step2 = 0;
            if type == 2
                for n = 1:size(obj,2)
                    QI = Q(step+1:step+obj(n).nCells);
                    p_1 = plot([step2 QI],[T1l(step+1) T1r(step+1:step+obj(n).nCells)]-273.15,'Color',[0 0.4470 0.7410]);
                    hold on
                    step = step + obj(n).nCells;
                    step2 = QI(end);
                end
            elseif type == 1
                p_1 = plot(QII,[T1l(1) T1r]-273.15,'Color',[0 0.4470 0.7410]);
                hold on
            end

            p_3 = plot(QII,[T3l(1) T3r]-273.15,'Color',[0.9290 0.6940 0.1250]);
            hold off
            xlabel("Q in MW")
            ylabel("T in °C")
            %FIX!
            % if obj(1).process_medium == "water" 
            %     legend(p_1 ,{"water", 'sand'})
            % elseif obj(1).process_medium == "air"
            %     legend(p_1,{"water", 'air'})
            % elseif obj(1).process_medium == "CO2"
            %     legend(p_1,{"water", 'CO2'})
            % end
            if obj(1).process_medium1 == "CO2"
                if obj(1).process_medium3 == "CO2"
                    legend({'CO_2', 'CO_2'})
                end
            end
            
        end
    end
    
    % density characteristic
    methods (Static)
        function rho=rho(p,h)
            spec = IF97;
            hsatl=spec.h(p,NaN,0);
            hsatv=spec.h(p,NaN,1);
            is2phase=hsatl<=h & h<=hsatv;
            
            x=NaN(size(h));
            x(is2phase)=(h(is2phase)-hsatl(is2phase))./(hsatv(is2phase)-hsatl(is2phase));
            
            T=NaN(size(h));
            T(~is2phase)=spec.T_ph(p(~is2phase),h(~is2phase));
            
            rho=1./spec.v(p,T,x);

        end

        function mf = mflu(dp,rhop,p,T,A,flu)
            uf = FluBed.wmf(dp,rhop,p,T);
            rho_air = DryAir.rho(p,T);
            mf = rho_air*uf*A*flu;
        end

    end
 end




