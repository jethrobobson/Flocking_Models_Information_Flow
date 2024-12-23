classdef Viscek_PP
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here

    properties

        base_dir
        save_dir

        eta
        N
        N1
        N2
        L
        dt
        r
        r11
        r12
        r21
        pos
        theta
        s
        v
        total_steps

        g
        omega
        leaders_inf


        delv_min
        delv_plus
        p_err
        nbins
        is
        zn
        zn_plus
        zn_neighb
        zn_neighb_plus
        select_freq
        count_
        elligible_zn_idx
        elligible_zn_check
        zn_plus2
        zn_neighbs2
        zn_plus2_neighb

        te
        leaders_te
        te_local
        te_peerage
        te_proletariat

        run_num

        holdoff_steps
        T_sim
        data_

        Cs
        Cc
        Cs_r
        Cc_r
    end

    methods
        function obj = Viscek_PP(eta, g, Npeerage, run_num)
            %Viscek_PP Construct an instance of this class
            %   INPUT:
            %       eta:    noise in Viscek Model [0,1]
            %       g:      true - peerage (big birds) are those with
            %       information on some resource in the distance
            %               false - peerage follows coupling dynamics
            %               (Explained below)
            %       Npeerage: Number of Peerage boids (Big Birds)
            %               if Npeerage = 0, reduces to standard Viscek
            %               model
            %       run_num:ensemble number for simulations. Also seeds the
            %               random number generator
            %   
            %% Hierarchical Vicsek Model
            % The flock has a two level heirarchy - the peerage and the proletariat.
            % Both can be thought of as sub-flocks. Members of the peerage pay
            % attention to one another with probability r11, irrespective of distance.
            % They also pay attention to the proletariat geometrically with probability
            % r12. The proletariat are a standard Vicsek model, with a geometric
            % connectivity. Furthermore the proletariat pay attention to the peerage
            % geometrically, as well as distance independently with probability r21.

            % Run this to see a Heirarchical Viscek Model in action.

            obj.base_dir = "P:/Research/Information_Flocking/Couzin/bin/Viscek_Peerage_Proletariat";
            %obj.base_dir = "/scratch/skelty/Information_Flocking/Couzin/bin/Viscek_Peerage_Proletariat";
            obj.save_dir = "";
            obj.run_num = run_num;
            rng(obj.run_num);

            obj.eta = eta;  % Individuality parameter set between 0 and 1
            obj.N = 200;    % Number of birds in the population
            %Npeerage = 25;
            if Npeerage == 0
                obj.N1 = [];
            else
                obj.N1 = 1:Npeerage; % Indices of the Peerage
            end
            %obj.N1 = 1:npeerage;  
            obj.N2 = (Npeerage+1):obj.N;  % Indices of the Proletariat
            obj.L = 15;     % Size of the world (effective interaction range, since r = 1)
            obj.dt = 0.1;   % Size of a single timestep
            obj.r = 1;  % Radius of Interaction

            if (g == false) && (Npeerage ~= 0)   
                %Coupling
                obj.r11 = .01; % Peerage-Peerage long-distance coupling
                obj.r12 = .5;  % Peerage-Proletariat short-distance coupling
                obj.r21 = .01; % Proletariat-Peerage long-distance coupling
                
            else
                %SVM
                obj.r11 = 0; % Peerage-Peerage long-distance coupling
                obj.r12 = 1;  % Peerage-Proletariat short-distance coupling
                obj.r21 = 0; % Proletariat-Peerage long-distance coupling
            end
            obj.s = 1;  % Speed

            %
            % g = false;
            obj.leaders_inf = obj.N1;
            if g == true
                obj.omega = .5;
            else
                obj.omega = 0;
                obj.leaders_inf = [];
            end
            obj.g = [0,1];
            
            obj.holdoff_steps = 500;
            obj.T_sim = 2000;
            
            obj = obj.reset_flock();

            obj = obj.initialize_anal();

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%% INITIALIZE & UPDATE PARAMETERS %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj = reset_flock(obj)
            obj.pos = rand([obj.N,2])*obj.L;
            obj.theta = rand([obj.N,1])*2*pi;
            obj.v = obj.s*[cos(obj.theta) sin(obj.theta)];

            obj.total_steps = 0;

            obj.delv_min = zeros(obj.N,2);
            obj.delv_plus = zeros(obj.N,2);
            %zn = del_x,del_y,del_z,del_vx,del_vy,del_vz between neighbors
            obj.zn = zeros(0,4);
            obj.zn_plus = zeros(0,4);
            obj.zn_neighb = [];
            obj.zn_neighb_plus = [];
            obj.zn_plus2_neighb = [];

            obj.select_freq = 10;
            obj.p_err = .23;
            obj.te_peerage = zeros(1);
        end

        function obj = initialize_anal(obj)
            obj.is(obj.total_steps + 1) = 0;
            obj.te(obj.total_steps + 1) = 0;
            obj.te_peerage(obj.total_steps + 1) = 0;
            obj.te_proletariat(obj.total_steps + 1) = 0;
            
        end
        %{

        function obj = get_measures(obj)

            if obj.total_steps + obj.holdoff_steps < 1
                obj.te_peerage(obj.total_steps + 1) = 0;
                obj.te_proletariat(obj.total_steps + 1) = 0;
            else
                obj.te_peerage(obj.total_steps + 1) = mean(obj.te_local(1:obj.N1));
                obj.te_proletariat(obj.total_steps + 1) = mean(obj.te_local(obj.N1+1:obj.N2));

            end
        end
        %}



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%% UPDATE FLOCK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        function obj = update_flock(obj)
            nns = rangesearch(obj.pos,obj.pos,obj.r,'Distance',@obj.periodic_dist);
            neighbs = obj.update_neighbors(nns);
            
            %obj = obj.get_zn_props(neighbs);

            obj = obj.update_velocity(neighbs);
            obj.pos = mod(obj.pos + obj.v*obj.dt,obj.L);
            obj.total_steps = obj.total_steps + 1;
            
            %{
            %if obj.total_steps < obj.sim_steps - obj.holdoff_steps
            if obj.total_steps + obj.holdoff_steps < 1
                obj = obj.initialize_anal();
            else
                %obj = obj.make_graph(repel_n,orient_n,attract_n);

                %obj.is(obj.total_steps + 1) = obj.information_storage();

                obj = obj.znplus2();
                [obj.leaders_te,obj.te(obj.total_steps + 1),obj.te_local]...
                    = obj.transfer_entropy2();

                obj = obj.get_measures();
                
            end
            %}
        end

        function obj = sim(obj)

            obj.save_dir = "N_peerage_"+string(length(obj.N1));
            if obj.omega == 0
                if isempty(obj.N1)
                    g = "/SVM";
                else
                    g = "/peerage_true";
                end
            else
                g = "/leaders_true";
            end
            file_ext = g + "/noise_"+obj.eta;
            obj.check_directory(file_ext);
            filename = obj.base_dir+ "/" + obj.save_dir + file_ext + "/data_"+obj.run_num + ".mat";
          
            tic;
            if ~exist(filename)
                obj = obj.reset_flock();
                if obj.holdoff_steps > 0
                    for i = (1:obj.holdoff_steps)
                        obj = obj.update_flock();
                    end
                end
    
                data(:,:,1) = [obj.pos angle(obj.v(:,1)+1i*obj.v(:,2))];
                %data(:,1) = [angle(obj.v(:,1)+1i*obj.v(:,2))];
                for tt = 2:obj.T_sim
                    obj = obj.update_flock();
                    %data(:,tt) = [angle(obj.v(:,1)+1i*obj.v(:,2))];
                    data(:,:,tt) = [obj.pos angle(obj.v(:,1)+1i*obj.v(:,2))];
                end
                obj.data_ = data;

                
            
                run_num = obj.run_num;
                eta = obj.eta;
                save(filename,"data","run_num","eta");
                %parsave(['Data/timeseries',num2str(ii,'%0.2u')],...
                %    data,ii)
                %disp([num2str(ii),' completed in ',num2str(toc),'s.'])
                obj.data_ = data;
                obj = obj.get_inf_measures(false);
            else
                disp("Simulation already run. Rerunning for non-symmetric case")
                f = load(filename); 
                obj.data_ = f.data;
                disp("Read in Data from: " +filename);
                obj.get_inf_measures(true);
            end
            toc;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function d = periodic_dist(obj,X,Y)
            %d = abs(X - Y);
            %d = abs(obj.L*(d > obj.L/2) - d);
            d = obj.periodic_rel_pos(X,Y);
            d = vecnorm(d,2,2);
        end

        function d = periodic_rel_pos(obj,X,Y)
            d = X - Y;
            d = d - obj.L*sign(d).*(abs(d) > obj.L/2);
        end

        function neighbs = update_neighbors(obj,nns)
            neighbs = {1,obj.N};
            for i = (1:obj.N)
                % PEERAGE (Big Birds)
                if i <= length(obj.N1)
                    % Neighbors are all big birds within the neighborhood,
                    % some big birds outside the neighborhood with prob r11
                    % , and some small birds within the neighborhood with
                    % prob r12
                    neighbs_ = union(cell2mat(nns(i)),find(rand([1,length(obj.N1)])<=obj.r11));
                    dummy = [];
                    for nn = neighbs_
                        if (nn > length(obj.N1))
                            if rand()<=obj.r12
                                dummy = [dummy nn];
                            end
                        else
                            dummy = [dummy nn];
                        end
                    end
                    %if isempty(dummy)
                    %    neighbs(i) = {[]};
                    %else
                        neighbs(i) = mat2cell(dummy,1);
                    %end

                % PROLETARIAT (Small Birds)
                else
                    dummy = union(cell2mat(nns(i)),find(rand([1,length(obj.N1)])<=obj.r21));
                    %if isempty(dummy)
                    %    neighbs(i) = {[]};
                    %else
                        neighbs(i) = mat2cell(dummy,1);
                    %end
                end
            end
        end

        function obj = update_velocity(obj,neighbs)
            new_v = zeros([obj.N,2]);
            %new_v = obj.v;
            for i = (1:obj.N)
                for nn = cell2mat(neighbs(i))
                    obj.count_ = obj.count_ + 1;
                    doop_s = obj.periodic_rel_pos(obj.pos(nn,:),obj.pos(i,:));
                    new_v(i,:) = new_v(i,:) + obj.v(nn,:);
                    if obj.elligible_zn_check(obj.count_) == 1
                        obj.zn_plus = cat(1,obj.zn_plus,cat(2,doop_s,obj.v(i,:) - obj.v(nn,:)));
                        obj.zn_neighb_plus = [obj.zn_neighb_plus i];
                        obj.zn_plus2_neighb = [obj.zn_plus2_neighb nn];
                    end
                end
                obj.delv_plus(i,:) = obj.v(i,:) - obj.s*new_v(i,:)/norm(new_v(i,:));

            end

            if obj.omega ~= 0
                new_v(obj.leaders_inf,:) = new_v(obj.leaders_inf,:)./vecnorm(new_v(obj.leaders_inf,:),2,2);
                new_v(obj.leaders_inf,:) = ((1-obj.omega)*new_v(obj.leaders_inf,:) + obj.omega*obj.g)/...
                        norm((1-obj.omega)*new_v(obj.leaders_inf,:) + obj.omega*obj.g);
            end

            theta = angle(new_v(:,1)+1i*new_v(:,2)) + pi*obj.eta*(2*rand([obj.N,1])-1);


            obj.v = obj.s*[cos(theta),sin(theta)];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%% TRANSFER ENTROPY %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj = check_config_ent_(obj,n_bins_)
        
            data = obj.data_;
            t = size(data,3);
            n = size(data,1);
            n1 = length(obj.N1);
            
        
            all_pos = data(1:n,1:2,1:t);
            pos_mean = mean(all_pos,1);
            all_pos = squeeze(all_pos - pos_mean);
        
            all_v = [cos(data(1:n,3,1:t)) sin(data(1:n,3,1:t))];
            %all_v = reshape(all_v,[n,t,2]);
            v_mean = mean(all_v,1);
            all_v = squeeze(all_v - v_mean);

            
        
        
            n_bins = n_bins_;
            r_bins = linspace(-obj.L/2,obj.L/2,n_bins);
            C = zeros(n_bins-1,n_bins-1,t,3);
            C_r = zeros(n_bins-1,n_bins-1,t,3);
            r_hist = zeros(n_bins-1,n_bins-1,t,3);
            
            for i = (1:t)
                for j = (1:2)
                    if j == 1
                        if n1 > 0
        
                            pos = squeeze(all_pos(1:n1,:,i));
                            v = squeeze(all_v(1:n1,:,i));
                            V = v*v';
                            P = pos*pos';
                        else
                            continue;
                        end
                    else
                        pos = squeeze(all_pos(n1+1:end,:,i));
                        v = squeeze(all_v(n1+1:end,:,i));
                        V = v*v';
                        P = pos*pos';
                    end
                        n = length(pos);
                    
                        %Get matrices of \Delta_x,\Delta_y, where each
                        %point in the grid is the x/y distance between
                        %bird i and bird j
                        
                        Rxy = reshape(kron(pos,ones([length(pos),1])),[n,n,2]);
                        RxyT = cat(3,Rxy(:,:,1)',Rxy(:,:,2)');
                        rxy = Rxy - RxyT;
                        rxy = rxy - obj.L*sign(rxy).*(abs(rxy) > obj.L/2);
                
                        % The bins to put correlations 
                        rxy = discretize(rxy,r_bins);
                    
                        for r_row = (1:n)
                            for r_col = (1:n)
                                C(rxy(r_row,r_col,1),rxy(r_row,r_col,2),i,j) = ...
                                    C(rxy(r_row,r_col,1),rxy(r_row,r_col,2),i,j)...
                                    + V(r_row,r_col);
                                C_r(rxy(r_row,r_col,1),rxy(r_row,r_col,2),i,j) = ...
                                    C_r(rxy(r_row,r_col,1),rxy(r_row,r_col,2),i,j)...
                                    + P(r_row,r_col);
                                r_hist(rxy(r_row,r_col,1),rxy(r_row,r_col,2),i,j) = ...
                                    r_hist(rxy(r_row,r_col,1),rxy(r_row,r_col,2),i,j)...
                                    + 1;
                            end
                        end 
                end
            end
            C(:,:,:,3) = sum(C,4); C_r(:,:,:,3) = sum(C_r,4);
            
            %r_hist(:,:,:,3) = sum(r_hist,4);
            %C = C./r_hist;
            C = squeeze(mean(C,3)); C_r = squeeze(mean(C_r,3));
            %r_hist = squeeze(mean(r_hist,3));
        
            
        
            f_k = zeros([n_bins-1,n_bins-1,3]);
            f_k_r = zeros([n_bins-1,n_bins-1,3]);
            for kx = (0:n_bins-2)
                for ky = (0:n_bins-2)
                    for x = (1:n_bins-1)
                        for y = (1:n_bins-1)
                            dum = exp((1i*2*pi*kx/(n_bins-1))*(x-(n_bins/2)))*...
                                exp((1i*2*pi*ky/(n_bins-1))*(y-(n_bins/2)));
                                for i = (1:3)
                                    f_k(kx+1,ky+1,i) = f_k(kx+1,ky+1,i) +dum*C(x,y,i);
                                    f_k_r(kx+1,ky+1,i) = f_k_r(kx+1,ky+1,i) +dum*C_r(x,y,i);
                                end
                        end
                    end
                end
            end
        
            spec = real(f_k);
            spec_r = real(f_k_r);
            %for i = (1:3)
            %    spec(:,:,i) = fftshift(spec(:,:,i));
            %end
            %spec(:,:,1)
            r_plot = (r_bins(:,1:end-1) + r_bins(:,2:end)) / 2;
        
        
            close all
            %{
            fig = figure('Position', [0 0 400 800]);
            ax1 = axes(fig,'Position',[0.1 .1 .8 .35]);
            ax2 = axes(fig,'Position',[0.1 0.55 .8 .35]);
            [X,Y] = meshgrid(r_plot,r_plot);
            %Z = rot90(fliplr(squeeze(IMs_avg(3,:,T_obs_id,4:end))),1);
            Z = spec;
            size(Z)
            size(X)
            size([X,Y]);
        
            mesh(ax1,X,Y,C(:,:,3))
            ylabel(ax1,'\Delta y') 
            xlabel(ax1,'\Delta x') 
            title(ax1,'Correlation Function')
        
        
            mesh(ax2,X,Y,r_hist(:,:,3))
            ylabel(ax2,'\Delta y') 
            xlabel(ax2,'\Delta x') 
            title(ax2,'Histogram of Distances')
        
        
            [X,Y] = meshgrid((0:n_bins-2),(0:n_bins-2));
            mesh(ax2,X,Y,spec(:,:,3))
            ylabel(ax2,'k_y') 
            xlabel(ax2,'k_x') 
            title(ax2,'Power Spectrum')
            zscale(ax2,'log')
            %}
            

            obj.Cs = zeros([3,1]);
            obj.Cc = zeros([3,1]);
            obj.Cs_r = zeros([3,1]);
            obj.Cc_r = zeros([3,1]);
            for i = 1:3
                modes_sqr = power(spec(:,:,i),2);
                pk_s = modes_sqr/sum(modes_sqr,'all');
                pk_c = modes_sqr/max(modes_sqr,[],'all');   

                obj.Cs(i) = -1*sum(pk_s.*log2(pk_s + eps),'all');
                obj.Cc(i) = -1*sum(pk_c.*log2(pk_c + eps),'all');


                modes_sqr = power(spec_r(:,:,i),2);
                pk_s = modes_sqr/sum(modes_sqr,'all');
                pk_c = modes_sqr/max(modes_sqr,[],'all');   

                obj.Cs_r(i) = -1*sum(pk_s.*log2(pk_s + eps),'all');
                obj.Cc_r(i) = -1*sum(pk_c.*log2(pk_c + eps),'all');
            end
            
        end

        function obj = get_inf_measures(obj,rerun)

        %% Information Measures
        % This section computes the mutual information measures as a function of
        % the observation length, number of binning states, and individuality
        % parameter. Bootstrapping the data made in Experiment 1, the code creates
        % an ensemble of observations, from each of which a distribution of mutual
        % informations is computed. 
        %% 1 Compute Mutual Information
        % Here we take the average angle made by the peerage and the proletariat
        % over experiment 1, X and Y respectively. Then we compute 4 different
        % mutual informations:
        %  1. How the immediate future of X depends on the present of X and Y
        %  2. How the immediate future of X depends on the present of X
        %  3. How the immediate future of Y depends on the present of X and Y
        %  4. How the immediate future of Y depends on the present of Y
        
        function AB = outer(A,B)
            AB = squeeze(reshape(A(:) * B(:).', [size(A), size(B)]));
        end

        
        
        
        %etas = linspace(0,1,52); etas = etas(2:(end-1));
        %Ts = kron([100 1000],1:9); Ts = [Ts 10000]; % Observation Lengths
        Ts = [10 kron([100],1:9) 1000];
        Ts = Ts(Ts < obj.T_sim); % Observation Lengths
        
        bin_width = [180 120 90 72 60 45 30 20 10 9 8 6 5 3 2];
        Ns = 360./bin_width; % Number of States
        Nens = 10;
        %base_dir = "P:/Research/Information_Flocking/Couzin/bin/Viscek_Peerage_Proletariat/";
        obj.save_dir = "/N_peerage_"+string(length(obj.N1));
        
        %filename = base_dir+ "/" + save_dir + "/"+obj.run_num + ".mat";
        filedir = obj.base_dir + obj.save_dir ;
        IMs = zeros([2,Nens,length(Ns),length(Ts)]);
        %obj.data_

        if length(obj.N1) > 0

            %Avg Direction per timestep of Peerage and Proletariat
            theta = [angle(squeeze(mean(exp(1i*obj.data_(obj.N1,3,:)),1)))...
                 angle(squeeze(mean(exp(1i*obj.data_(obj.N2,3,:)),1)))]';
            
        
                theta = theta + pi;
            
            %loop over the state coarsness
                for nn = 1:length(Ns)
                    % bin the timeseries
                    stateTimeSeries = discretize(theta,linspace(-pi,pi,Ns(nn)+1));

                    %disp([obj.eta nn/length(Ns) toc/60])
                    % loop over the observation times
                    for tt = 1:length(Ts)
                        Tobs = Ts(tt);
                        
                        % bootstrap to create an ensemble
                        for jj = 1:Nens
                            %disp([ii nn/length(Ns) tt/length(Ts) jj/Nens toc/60])
                            % pick a random subsequence from the data
                            idx = (1:Tobs) + randi(obj.T_sim-Tobs);
                            X = stateTimeSeries(1,idx);
                            Y = stateTimeSeries(2,idx);
                            % Populate the joint distribution with the subdata
                            pXXY = zeros(Ns(nn)*[1 1 1]);
                            pYYX = zeros(Ns(nn)*[1 1 1]);
                            for kk = 2:Tobs
                                pXXY(X(kk),X(kk-1),Y(kk-1)) = ...
                                    pXXY(X(kk),X(kk-1),Y(kk-1)) + 1;
                                pYYX(Y(kk),Y(kk-1),X(kk-1)) = ...
                                    pYYX(Y(kk),Y(kk-1),X(kk-1)) + 1;
                            end
                            % Incorporate rotational symmetry into distribution
                            temp1 = pXXY;
                            temp2 = pYYX;
                            if rerun == false
                                for kk = 1:(Ns(nn)-1)
                                    temp1 = temp1 + circshift(pXXY,kk*[1 1 1]);
                                    temp2 = temp2 + circshift(pYYX,kk*[1 1 1]);
                                end
                            end
                            % Normalize the distributions
                            pXXY = temp1/sum(temp1(:));
                            pYYX = temp2/sum(temp2(:));

                            % Get the marginal distributions
                            pXX = squeeze(sum(pXXY,3));
                            pXY = squeeze(sum(pXXY,2));
                            p_XY = squeeze(sum(pXXY,1));
                            pX_ = squeeze(sum(pXX,2));
                            p_X = squeeze(sum(pXX,1));
                            pYY = squeeze(sum(pYYX,3));
                            pYX = squeeze(sum(pYYX,2));
                            p_YX = squeeze(sum(pYYX,1));
                            pY_ = squeeze(sum(pYY,2));
                            p_Y = squeeze(sum(pYY,1));
                            

                            % M(X(t+dt):X(t)&Y(t))
                            pXpXY = outer(pX_,p_XY);
                            idx = pXXY.*pXpXY ~= 0;
                            temp = pXXY.*log2(pXXY./pXpXY);
                            IMs(1,jj,nn,tt) = sum(temp(idx),'omitnan');

                            % M(X(t+dt):X(t))
                            pXpX = outer(pX_,p_X);
                            idx = pXX.*pXpX ~= 0;
                            temp = pXX.*log2(pXX./pXpX);
                            IMs(2,jj,nn,tt) = sum(temp(idx),'omitnan');

                            % T(Y->X) = M(X(t+dt):X(t)&Y(t)) - M(X(t+dt):X(t))
                            % SMALL TO BIG
                            IMs(3,jj,nn,tt) = IMs(1,jj,nn,tt) - IMs(2,jj,nn,tt);
                            
                            % M(Y(t+dt):X(t)&Y(t))
                            pYpYX = outer(pY_,p_YX);
                            idx = pYYX.*pYpYX ~= 0;
                            temp = pYYX.*log2(pYYX./pYpYX);
                            IMs(4,jj,nn,tt) = sum(temp(idx),'omitnan');

                            % M(Y(t+dt):Y(t))
                            pYpY = outer(pY_,p_Y);
                            idx = pYY.*pYpY ~= 0;
                            temp = pYY.*log2(pYY./pYpY);
                            IMs(5,jj,nn,tt) = sum(temp(idx),'omitnan');

                            % T(X->Y) = M(Y(t+dt):X(t)&Y(t)) - M(Y(t+dt):Y(t))
                            % BIG TO SMALL
                            IMs(6,jj,nn,tt) = IMs(4,jj,nn,tt) - IMs(5,jj,nn,tt);

                            %{
                            % T(Y->X) = H(X|X) - H(X|X'Y) (Same as I3)
                            %H(X|X') = sum(p(XX')log(p(X,X')/p(X'))
                            logdenom = outer(ones([1,Ns(nn)]),p_X);
                            idx = pXX.*logdenom ~= 0;
                            temp = pXX.*log2(pXX./logdenom);
                            HXX = -sum(temp(idx),'omitnan');
                            %HX = -sum(pX_.*log2(pX_),'omitnan');

                            %H(X|X'Y) = sum(p(XX'Y)log(p(XX'Y)/p(X'Y))
                            logdenom = outer(ones([1,Ns(nn)]),p_XY);
                            idx = pXXY.*logdenom ~= 0;
                            temp2 = pXXY.*log2(pXXY./logdenom);
                            HXXY = -sum(temp2(idx),'omitnan');
                            
                            IMs(7,jj,nn,tt) = HXX - HXXY;

                            % T(X->Y) = H(Y|Y) - H(Y|Y'X) (same as I6)                          
                            %H(Y|Y') = sum(p(YY')log(p(YY')/p(Y'))
                            logdenom = outer(ones([1,Ns(nn)]),p_Y);
                            idx = pYY.*logdenom ~= 0;
                            temp = pYY.*log2(pYY./logdenom);
                            HYY = -sum(temp(idx),'omitnan'); 
                            %HY = -sum(pY_.*log2(pY_),'omitnan')
                            
                            %H(Y|Y'X) = sum(p(YY'X)log(p(YY'X)/p(Y'X))
                            logdenom = outer(ones([Ns(nn),1]),p_YX)
                            idx = pYYX.*logdenom ~= 0;
                            pYYX./logdenom
                            temp2 = pYYX.*log2(pYYX./logdenom);
                            HYYX = -sum(temp2(idx),'omitnan');
                            IMs(8,jj,nn,tt) = HYY - HYYX

                            %HY - HYYX %Same as I4
                            %HY - HYY %Same as I5
                            %}

                            
                        end
                    end
                end

        end

        n_bins_ = 20;
        obj = obj.check_config_ent_(n_bins_);
            1
            flock_dir_t = mean(exp(1i*obj.data_(:,3,:)),1);
            theta_avg = mean(abs(flock_dir_t));
            flock_speed_sqr_avg = mean(abs(flock_dir_t).^2);
            delta_theta = flock_speed_sqr_avg - theta_avg^2;

            big_flock_dir_t = mean(exp(1i*obj.data_(obj.N1,3,:)),1);
            theta_big_avg = mean(abs(big_flock_dir_t));
            flock_speed_sqr_avg = mean(abs(big_flock_dir_t).^2);
            delta_theta_big = flock_speed_sqr_avg - theta_big_avg^2;

            small_flock_dir_t = mean(exp(1i*obj.data_(obj.N2,3,:)),1);
            theta_small_avg = mean(abs(small_flock_dir_t));
            flock_speed_sqr_avg = mean(abs(small_flock_dir_t).^2);
            delta_theta_small = flock_speed_sqr_avg - theta_small_avg^2;

            run_num = obj.run_num; eta = obj.eta;
            config_ent = obj.Cs; config_comp = obj.Cc;
            config_ent_r = obj.Cs_r; config_comp_r = obj.Cc_r;
            
            if obj.omega == 0
                if isempty(obj.N1)
                    g = "/SVM";
                else
                    g = "/peerage_true";
                end
            else
                g = "/leaders_true";
            end
            file_ext = g + "/noise_"+obj.eta;
            obj = obj.check_directory(file_ext);
            if rerun == true
                i = "/IM_no_sym_";
            else
                i = "/IM_";
            end
            filename = filedir + file_ext+i+obj.run_num+".mat";
            save(filename,"run_num","eta","IMs","Nens","Ns","Ts",...
                "theta_avg","delta_theta","theta_big_avg",...
                "delta_theta_big","theta_small_avg","delta_theta_small", ...
                "config_ent","config_comp","config_ent_r","config_comp_r");
            %parsave(['Data/infomeasures_',...
            %    num2str(ii,'%0.2u')],...
            %    IMs,ii);
            toc
        
        end

        function obj = rerun_inf_measures(obj)

            obj.save_dir = "N_peerage_"+string(length(obj.N1));
            if obj.omega == 0
                if isempty(obj.N1)
                    g = "/SVM";
                else
                    g = "/peerage_true";
                end
            else
                g = "/leaders_true";
            end
            file_ext = g + "/noise_"+obj.eta;
            obj.check_directory(file_ext);
            filename = obj.base_dir+ "/" + obj.save_dir + file_ext + "/data_"+obj.run_num + ".mat";
          
            tic;
            if exist(filename)
                f = load(filename); 
                obj.data_ = f.data;
                disp("Read in Data from: " +filename);
                obj.get_inf_measures(true);
            else
                disp(filename +": Does not exist. Run Simulation and Save Data.")
            end
        end


        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        function obj = update_all_and_plot(obj,num_steps)
            close all
            %fig = figure('Position',[0 0 600 600]);

            fig = figure('Position', [0 0 600 600]);
            %y_spac = .1;
            %quiv_pos = .3;
            ax1 = axes(fig,'Position',[0.1 .1 0.8 .8]);
            %ax1 = axes(fig,'Position',[0.1 (quiv_pos) 0.8 (1-quiv_pos-y_spac)]);
            %ax2 = axes(fig,'Position',[0.1 0.1 .8 y_spac]);

            %ax  = axes(fig,'Position',[0 0 1 1]);
            scaleFactor = 0.5;
            Peerage = quiver(ax1,...
                obj.pos(obj.N1,1)-scaleFactor*obj.v(obj.N1,1)*0.5,...
                obj.pos(obj.N1,2)-scaleFactor*obj.v(obj.N1,2)*0.5,...
                scaleFactor*obj.v(obj.N1,1),scaleFactor*obj.v(obj.N1,2),...
                'Color',[0.5 0 0.5],...
                'AutoScale','off',...
                'AutoScaleFactor',5/obj.L,...
                'LineWidth',2);
            hold(ax1,'on')
            Proletariat = quiver(ax1,obj.pos(obj.N2,1),obj.pos(obj.N2,2),...
                0.5*scaleFactor*obj.v(obj.N2,1),0.5*scaleFactor*obj.v(obj.N2,2),...
                'Color',[0.8 0 0],...
                'AutoScale','off',...
                'AutoScaleFactor',5/obj.L,...
                'LineWidth',1);
            set(ax1,...
                'XLim',[0 obj.L],...
                'YLim',[0,obj.L],...
                'XTick',[],...
                "YTick",[])
            ButtonHandle = uicontrol('Style', 'PushButton', ...
                'String', 'Exit', ...
                'Callback', 'delete(gcbf)');
            %title(ax1, "N = "+string(obj.N) + ", N_{big} = " + string(length(obj.N1))...
             %   + ", Noise \eta = " + string(obj.eta))
            title(ax1, "N = "+string(obj.N) + ", L = " + string(obj.L)...
                + ", Noise \eta = " + string(obj.eta*2*pi))
            hold(ax1,"off");
            %te_pee = obj.te_peerage(1:obj.total_steps+1);
            %te_pro = obj.te_proletariat(1:obj.total_steps+1);

            %x_  = linspace(0,obj.total_steps,obj.total_steps+1);
            %inf_ = plot(ax2,x_,te_pee,x_,te_pro);
            %legend(inf_, 'Peerage','Proletariat','Location','northeast')
            %ax1.Title.String = 'Num Resource Leaders = ' + string(length(obj.leaders_inf)) + ': N = '+string(obj.N)...
            %        +', Density = '+string(obj.N/(obj.L^2))+', \omega = ' + string(obj.omega)...
            %        +', \eta = ' + obj.eta;

            for i = obj.total_steps+1:(num_steps + obj.total_steps)
                obj = obj.update_flock();

                Peerage.XData = obj.pos(obj.N1,1)-scaleFactor*obj.v(obj.N1,1)*0.5;
                Peerage.YData = obj.pos(obj.N1,2)-scaleFactor*obj.v(obj.N1,2)*0.5;
                Peerage.UData = scaleFactor*obj.v(obj.N1,1);
                Peerage.VData = scaleFactor*obj.v(obj.N1,2);
                Proletariat.XData = obj.pos(obj.N2,1)-0.5*scaleFactor*obj.v(obj.N2,1)*0.5;
                Proletariat.YData = obj.pos(obj.N2,2)-0.5*scaleFactor*obj.v(obj.N2,2)*0.5;
                Proletariat.UData = 0.5*scaleFactor*obj.v(obj.N2,1);
                Proletariat.VData = 0.5*scaleFactor*obj.v(obj.N2,2);

                %x_ = linspace(0,obj.total_steps,obj.total_steps+1);
                
                %set(inf_(1),'XData',x_,'YData',obj.te_peerage(1:obj.total_steps+1));
                %set(inf_(2),'XData',x_,'YData',obj.te_proletariat(1:obj.total_steps+1));

                F(i) = getframe(fig) ;
                drawnow
                if ~ishandle(ButtonHandle)
                    disp('Loop stopped by user');
                    break;
                end
            end
            writerObj = VideoWriter(obj.base_dir + '/TE_pre_pro_test.avi');
            writerObj.FrameRate = 24;

            open(writerObj);

            for i=1:length(F)
                % convert the image to a frame
                frame = F(i) ;    
                writeVideo(writerObj, frame);
            end
            % close the writer object
            close(writerObj);
        end

        function obj = check_directory(obj,file_ext)
            basedir = obj.base_dir+"/"+obj.save_dir;

			dir1 = basedir+"/"+file_ext;

			%dir2 = basedir+"/del_ro"+string(obj.del_ro)+"/del_ra"+string(obj.del_ra);

            %disp(dir1);

            if ~exist(basedir, 'dir')
               mkdir(basedir)
            end
            if ~exist(dir1, 'dir')
               mkdir(dir1)
            end
            %if ~exist(dir2, 'dir')
            %   mkdir(dir2)
            %end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%% LIZIER TRANSFER ENTROPY %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj = get_zn_props(obj,neighbs)
            sum_neighbors_pee = 0;
            sum_neighbors_pro = 0;
            for id_ = (1:length(obj.N1))
                sum_neighbors_pee = sum_neighbors_pee + cellfun(@length,neighbs(id_));
            end
            for id_ = (length(obj.N1)+1:obj.N)
                sum_neighbors_pro = sum_neighbors_pro + cellfun(@length,neighbs(id_));
            end
            
            d1 = randsample((1:sum_neighbors_pee),min(max(obj.N,int16(sum_neighbors_pee/obj.select_freq)),sum_neighbors_pee),false);
            d2 = randsample((sum_neighbors_pee+1:sum_neighbors_pro+sum_neighbors_pee),...
                min(max(obj.N,int16(sum_neighbors_pro/obj.select_freq)),sum_neighbors_pro),false);

            obj.elligible_zn_idx = cat(2,d1,d2);
            
            obj.elligible_zn_check = zeros(sum_neighbors_pee + sum_neighbors_pro,1);
            obj.elligible_zn_check(obj.elligible_zn_idx) = 1;
            
            obj.elligible_zn_check = zeros(sum_neighbors_pee + sum_neighbors_pro,1);
            obj.count_ = 0;
            if obj.total_steps > 0
                obj.delv_min = obj.delv_plus;
                obj.zn = obj.zn_plus;
                obj.zn_plus = zeros(0,4);

                obj.zn_neighb= obj.zn_neighb_plus;
                obj.zn_neighb_plus = [];
                obj.zn_neighbs2 = obj.zn_plus2_neighb;
                obj.zn_plus2_neighb = [];
            end

        end


        function p = p_args(obj,args_,params_)
            % args_ = function values at which to compute probability
            % params_ = values describing state of the system
            % eg: for p_qn, params_ =  obj.delv_min
            %           args_ =  some choice of change in v, usually
            %           obj.delv_min(id,:)
            p = 0;
            for id = (1:size(params_,1))
                p = p + exp(obj.get_exp_arg_(args_,params_(id,:)));
            end
        end


        function exp_arg_ = get_exp_arg_(obj,args_,params_)
            exp_arg_ = 0;
            for i = (1:length(args_))
                exp_arg_ = exp_arg_ - ((args_(i) - params_(i))/obj.p_err)^2;
            end
        end

        function is = information_storage(obj)

            is = 0;
            %norm_ = 1;
            norm_ = obj.N;
            
            for q = (1:obj.N)
                pn = obj.p_args(obj.delv_min(q,:),obj.delv_min);
                pnplus = obj.p_args(obj.delv_plus(q,:),obj.delv_plus);
                pnnplus = obj.p_args(cat(2,obj.delv_min(q,:),obj.delv_plus(q,:)),...
                    cat(2,obj.delv_min,obj.delv_plus));

                is = is + log2(norm_ * (pnnplus+eps)/(pn*pnplus+eps));
            end
            
        end


        function obj = znplus2(obj)

            k = 1;
            obj.zn_plus2 = zeros(size(obj.zn_neighb,2),4);
            %size(obj.zn_neighb,1)
            %size(obj.zn_plus2_neighb,2)
            for i = obj.zn_neighb
                nn = obj.zn_neighbs2(k);
                del_s = obj.periodic_rel_pos(obj.pos(i,:),obj.pos(nn,:));
                obj.zn_plus2(k,:) = cat(2,del_s,obj.v(i,:) - obj.v(nn,:));
                k = k+1;
            end
        end
        
        function params_ = get_qnznznplus_params_(obj)
            params_ = [];
            k = 1;
            
            for i = obj.zn_neighb
                params_ = cat(1,params_,cat(2,obj.delv_min(i,:),obj.zn(k,:),obj.zn_plus2(k,3:4)));
                
                k = k+1;
            end
        end
        function params_ = get_znznplus_params_(obj)
            params_ = [];
            k = 1;
            
            for i = obj.zn_neighb
                params_ = cat(1,params_,cat(2,obj.zn(k,3:4),obj.zn_plus2(k,3:4)));
                k = k+1;
            end
        end


        function params_ = get_qnqnpluszn_params_(obj)
            params_ = [];
            k = 1;
            for i = obj.zn_neighb
                params_ = cat(1,params_,cat(2,obj.delv_min(i,:),obj.delv_plus(i,:),obj.zn(k,:)));
                k = k+1;
            end
        end

        function params_ = get_qnzn_params_(obj)
            params_ = [];
            k = 1;
            for i = obj.zn_neighb
                params_ = cat(1,params_,cat(2,obj.delv_min(i,:),obj.zn(k,:)));
                k = k+1;
            end
        end

        function params_ = get_qnzn2_params_(obj)
            params_ = [];
            k = 1;
            for i = obj.zn_neighb
                params_ = cat(1,params_,cat(2,obj.delv_min(i,:),obj.zn(k,:)));
                k = k+1;
            end
        end
        

        function [leaders,te,te_lst] = transfer_entropy2(obj)
            %bins = linspace(-obj.s,obj.s,6);
            %te_lst
            te_lst = zeros(obj.N,1);
            
            qnznznplus_params = obj.get_qnznznplus_params_();
            znznplus_params = obj.get_znznplus_params_();
            qnzn_params = obj.get_qnzn2_params_();
            %norm_ = (pi()* obj.p_err^2)^3/2;
            norm_ = 1;
            for q = (1:obj.N)
              
                pqnznznplus = 0;
                pznznplus = 0;
                pqnzn = 0;
                pzn = 0;

                zn_ind = find(obj.zn_neighb == q);

                for i = zn_ind
                    pzn = pzn + obj.p_args(obj.zn(i,3:4),obj.zn(:,3:4));

                    args_ = cat(2,obj.delv_min(q,:),obj.zn(i,:));
                    pqnzn = pqnzn + obj.p_args(args_,qnzn_params);

                    args_ = cat(2,obj.zn(i,3:4),obj.zn_plus2(i,3:4));
                    pznznplus = pznznplus + obj.p_args(args_,znznplus_params);
                    
                    args_ = cat(2,obj.delv_min(q,:),obj.zn(i,:),obj.zn_plus2(i,3:4));
                    pqnznznplus = pqnznznplus + obj.p_args(args_,qnznznplus_params);
                end

                te_lst(q) = log2(norm_*(pqnznznplus * pzn  + eps)/(pznznplus*pqnzn + eps));
            end
            te = sum(te_lst);
            [te_,leaders] = maxk(te_lst,10);

        end

        
    end
    
end