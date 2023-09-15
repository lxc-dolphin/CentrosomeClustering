

# Packages are using
using LinearAlgebra
using DifferentialEquations
using Distances
using Combinatorics
using Plots
using Clustering
using Interpolations


# Radial force and energy
"""
X: target CT position to simulate, force is applied on X
Y: position on the geo which is the shortest distance to X
dis: distance between X and Y
K: intensity of radial force
P: power to distance dis^P, default: 1
minL: minimum length/threshold, default: 0
"""

function RadialForce(dis,X,Y,K,P) # direction is already applied
    # r = maximum([dis,minL])
    r = dis
    if norm(Y-X) != 0
        force_ra = -K*r^P*(X-Y)/norm(X-Y)
    else
        force_ra = [0;0;0]
    end

    return force_ra
end

function RadialEnergy(dis,K)

    r = dis

    RaEnergy = 0.5*K*r^2

    return RaEnergy
end

# Shortest distance from a point to ellipsoid surface
"""

The shortest distance from a point to Triaxial Ellipsoid or Biaxial Ellipsoid or Sphere

  (x/a)^2+(y/b)^2+(z/c)^2=1   Triaxial Ellipsoid Equation centered at the
   origin
 Parameters:
 * X, [x y z]     - Cartesian coordinates data, n x 3 matrix Lxc!199011or three n x 1 vectors
 * axis,[a; b; c] - ellipsoid radii  [a; b; c],its axes % along [x y z] axes

                  For Triaxial ellipsoid ,it must be a > b > c

                  For Biaxial ellipsoid ,it must be a = b > c

                 For Sphere ,it must be a = b = c

 Output:
 * Xo,[xo yo zo]   -  Cartesian coordinates of Point onto ellipsoid
"""
function ProjDisGeo(X,axis)

eps = 0.0001
a, b, c  = axis

x, y, z = X
lambda = 0.1

E, F, G  = sign(a)/a^2, sign(b)/b^2, sign(c)/c^2

xo, yo, zo = a*x/sqrt(x^2+y^2+z^2), b*y/sqrt(x^2+y^2+z^2), c*z/sqrt(x^2+y^2+z^2)

if abs(x^2/a^2+y^2/b^2+z^2/c^2-1)<=eps

    xo = x;
    yo = y;
    zo = z;

else

    for i=1:50
        j11=1+lambda*E;
        j14=xo*E;
        j22=1+lambda*F;
        j24=yo*F;
        j33=1+lambda*G;
        j34=zo*G;
        j41=2*xo*E;
        j42=2*yo*F;
        j43=2*zo*G;
        A=[ j11   0    0   j14;
             0   j22   0   j24;
             0    0   j33  j34;
            j41  j42  j43   0];

        f1=(xo-x)+xo*lambda*E;
        f2=(yo-y)+yo*lambda*F;
        f3=(zo-z)+zo*lambda*G;
        f4=E*xo^2+F*yo^2+G*zo^2-1;
        de_f = [f1 f2 f3 f4]';


        if maximum(abs.(de_f))>=eps
            de_x=-A\de_f;
        else
            de_x = zeros(4,1);
        end
        xo=xo+de_x[1];
        yo=yo+de_x[2];
        zo=zo+de_x[3];
        lambda=lambda+de_x[4];

    end

end

Y = [xo,yo,zo]
dis=sqrt((x-xo)^2+(y-yo)^2+(z-zo)^2)
    return Y, dis

end

# intercentrosomal force
"""
Calculate intercentrosomal force
dis: distance between X and Y
X: target CT position, force will be applied on X
Y: other CT position
K1,K2: attraction, repulsion intensity
L1,L2: effective widthes of attraction and repulsion
PL: effective range of InterCT
"""
function InterCTForce(dis,K1,K2,L1,L2,PL) # return scalar
    if dis >= PL
        mag = 0.0
    elseif dis <= 0.05 # um
        mag = 0.0
    else
        mag = -(K1*exp(-dis/L1)-K2*exp(-(PL-dis)/L2))
    end

    return mag
end


function InterCTenergy(dis,K1,K2,L1,L2,PL)
    if dis >= PL
        CTenergy = K1*L1*(1-exp(-PL/L1))
    elseif dis <= 0.05 # um
        CTenergy = 0.0
    else
        CTenergy = K1*L1*(1-exp(-dis/L1)) + K2*L2*(1-exp(-(PL-dis)/L2))
    end

    return CTenergy
end



function CT_cluster(du,u,p,t)

    k_a, k_r, l_a, l_r, p_l, k_in, k_out, a, b, c, gamma, xi, D, num_CT = p

    # du = []*length(u) # initial du


    ind_CT = [Int(i) for i in 1:num_CT]
    posi_CT =[u[(Int(i)-1)*3+1:(Int(i)-1)*3+3] for i in 1:num_CT]
    dis_mat = pairwise(Euclidean(),posi_CT)
    F_c_mat = InterCTForce.(dis_mat,k_a, k_r, l_a, l_r, p_l)

    # choose which CT to pair with the target CT, cuurently not used
    pairs = combinations(ind_CT,Int(num_CT-1))|> collect
    # choose which two CTs are paired
    CT_pairs = combinations(ind_CT,2)|> collect

    # inter-CT force fluctuation
    for i in 1:length(CT_pairs)

        m, n = CT_pairs[i]
        F_c = F_c_mat[m,n]
        du[Int(i+num_CT*3)] = -xi*(u[Int(i+num_CT*3)]-F_c)

    end

    # CT position
    for u_tar in posi_CT

        du_tar = zeros(3) # vector 3*1; initial du_tar
        ind_tar = findall(e->e==u_tar,posi_CT)
        #posi_CT_other = posi_CT
        posi_CT_other = [e for e in posi_CT] # left other CTs
        deleteat!(posi_CT_other, ind_tar)

        # interCT force
        for u_other in posi_CT_other


            ind_other = findall(e->e==u_other,posi_CT)
            dis = dis_mat[ind_tar,ind_other]
            F_c = F_c_mat[ind_tar,ind_other]

            ind_fcf = findall(e->e == sort([ind_tar;ind_other]), CT_pairs)
            F_cf = u[Int(ind_fcf[1]+num_CT*3)]


            du_tar = du_tar + F_cf/gamma*(u_tar-u_other)/dis

        end
        # radial force
        (x_o, dis_o) = ProjDisGeo(u_tar,[a,b,c])

        if  ( (u_tar[1]/a)^2 + (u_tar[2]/b)^2 + (u_tar[3]/c)^2 )> 1
            F_ra = RadialForce(dis_o,u_tar,x_o,k_in,1)
        else
            F_ra = RadialForce(dis_o,u_tar,x_o,k_out,1)
        end
        du_tar = du_tar + F_ra/gamma

        du[(ind_tar[1]-1)*3+1:(ind_tar[1]-1)*3+3] = du_tar
    end

end


# noise term

dσ = function CT_cluster_noise(dσ,u,p,t)
    KbT_ = 4.1E-6 # nN*um
    k_a, k_r, l_a, l_r, p_l, k_in, k_out, a, b, c, gamma, xi, D, num_CT, x_eta = p


    eta = 2*KbT_/gamma

    eta = (x_eta^2)*eta

    dσ[1:Int(num_CT*3)] = sqrt(eta)*ones(Int(num_CT*3))
    dσ[1+Int(num_CT*3):end] = sqrt(2*D)*ones(Int(num_CT*(num_CT-1)/2))

end
# Randomly initialize centrosome position
function InitCTPosi(p)
    k_a, k_r, l_a, l_r, p_l, k_in, k_out, a, b, c, gamma, xi, D, num_CT = p
    #  a,b,c and num_CT will be used
    rela_dis = 0
    ind_CT = [Int(i) for i in 1:num_CT]
    CT_pairs = combinations(ind_CT,2)|> collect
    dis_seri_out = Float64[]
    X_init_out = Float64[]

    LowBound = 5
    if num_CT>4
        LowBound = 2.5
    end

    while (rela_dis < LowBound)
        num_CT = Int(num_CT)
        x = rand(Float64, (num_CT,1))
        phi = x*2*pi
        y = rand(Float64, (num_CT,1))
        theta = y*pi

        X_init = Float64[]

        for i in 1:num_CT
            x_temp = para_xyz(a,b,c,theta[i],phi[i])

            X_init = [X_init;x_temp]
        end
        X_init_out  = X_init
        posi_CT_temp =[ X_init[(Int(i)-1)*3+1:(Int(i)-1)*3+3] for i in 1:num_CT]
        dis_mat = pairwise(Euclidean(),posi_CT_temp)

        # dis_seri_out = dis_mat

        dis_seri_out = [ dis_mat[n,m] for (n,m) in CT_pairs]
        rela_dis = minimum(dis_seri_out)
    end

    F_c_seri = InterCTForce.(dis_seri_out,k_a, k_r, l_a, l_r, p_l)
    # F_c_seri = [ F_c_mat[n,m] for (n,m) in CT_pairs]

    return X_init_out, F_c_seri

end

"""
Functions for the clustering algrithm

eg. clustering number/pole; clustering time; orientation....
"""
# orientation, calculate the angle between the pole vector and x-y plane
function Orien_angle(Posi, p)

    ang_temp = 0
    orient_seri = Vector{Float64}()

    k_a, k_r, l_a, l_r, p_l, k_in, k_out, a, b, c, gamma, xi, D, num_CT = p
    for i = 1:Int(num_CT)
        po_each = Posi[i] # x, y, z locations
        po_vect = po_each/norm(po_each)

        ang_temp = ang_temp + asin(po_vect[3]/1)/(pi/2)*90
        append!(orient_seri,asin(po_vect[3]/1)/(pi/2)*90)
    end

    Mean_anlge = ang_temp/num_CT

    return Mean_anlge, orient_seri
end




# clasification for each step
function NaiveClas(posi)
    dis_0 = pairwise(Euclidean(),posi)


    link0 = hclust(dis_0,linkage =:single)
    cut0 = cutree(link0, h=3) # 3um threshold

    return maximum(cut0), cut0
end


# get the position seri data from solutions
function Get_Posi_seri(sol_u,num_CT)
    totT = size(sol_u)[1]

    Posi_seri = []
    for i in 1:totT
        posi = sol_u[i][1:3*num_CT]
        posi_CT_each =[posi[(Int(i)-1)*3+1:(Int(i)-1)*3+3] for i in 1:num_CT]
        Posi_seri = [Posi_seri; [posi_CT_each]]
    end

    return Posi_seri
end
# get the position X Y X seris from solutions
function Get_XYZ_seri(sol_u,num_CT)
    totT = size(sol_u)[1]
    X_seri = zeros(num_CT,2);Y_seri = zeros(num_CT,2);Z_seri = zeros(num_CT,2)


    for i in 1:totT
        x_each = zeros(num_CT);y_each = zeros(num_CT);z_each = zeros(num_CT)
        for j in 1:num_CT
            x_each[j] = sol_u[i][(j-1)*3+1]
            y_each[j] = sol_u[i][(j-1)*3+2]
            z_each[j] = sol_u[i][(j-1)*3+3]
        end
        X_seri = [X_seri x_each]
        Y_seri = [Y_seri y_each]
        Z_seri = [Z_seri z_each]
    end

    return [X_seri[:,3:end], Y_seri[:,3:end], Z_seri[:,3:end]]
end


# clasification with position series
function NaiveClas_seri(Posi_seri)
    totT = size(Posi_seri)[1]

    N_clusTotal = zeros(totT) # vector
    Indx_clus = zeros(totT)
    for i in 1:totT
        N_clusTotal[i],Indx_clus[i] = NaiveClas(Posi_seri[i])
    end

    return N_clusTotal
end

function NaiveClas_seri1(Posi_seri)
    totT = size(Posi_seri)[1]
    CT_num = size(Posi_seri[1])[1]

    N_clusTotal = zeros(totT) # vector

    Indx_clus = zeros(CT_num,totT)
    for i in 1:totT
        N_clusTotal[i],Indx_clus[:,i] = NaiveClas(Posi_seri[i])
    end

    return N_clusTotal, Indx_clus
end
# Parameterize CT position from angles
function para_xyz(a,b,c,theta,phi)

    x = a*cos(theta)*cos(phi)
    y = b*cos(theta)*sin(phi)
    z = c*sin(theta)

    return [x,y,z]
end

# functions for 2 CTs and mulit CTs
"""
Calculate distance between 2 CTs
input: X,Y,Z ->
positions of 2 CTs: size(X) = [2,length of steps]
output: Dis ->
distance between 2CTs: size(Dis) = length of steps
"""

function Dis_2CT(X,Y,Z)
    Dis = sqrt.((X[1,:]-X[2,:]).^2
              + (Y[1,:]-Y[2,:]).^2
              + (Z[1,:]-Z[2,:]).^2)
    return Dis
end

"""
Cluster Time with 4 or Multi CTs
"""
function Clus_T(N_clusT,time)

    Thre_tf = 30 # threshhold time to filter any noisy mulitipolar
    # unit: second
    Thre_tb = 300 # threshold time to verify a acceptable bipolar state

    trusted_b = 0
    ClustT = 0

    interp_linear = LinearInterpolation(time, N_clusT)
    time_intp = 0:1:last(time)
    N_clus_intp = interp_linear(time_intp)

    # display(plot(time_intp,N_clus_intp))

    # find the indx of bipolar state,
    # the indx number is the time mark with unit second
    Indx_bi = findall(x -> x == 2, N_clus_intp)
    # count the time difference
    diff_indx_b = diff(Indx_bi)
    # if the difference is 1 which means the bipolar state is continuous
    #                larger than 1 means there are other states,
    #                filter out these states if it is no longer than Thre_tf
    #                in other words, larger than Thre_tf can be considered as other states



    overTh_index = findall(x -> x > Thre_tf, diff_indx_b)



    Indx_Not_bi = findall(x -> x != 2, N_clus_intp)
    diff_indx_nb = diff(Indx_Not_bi)
    overTh_index_nb = findall(x -> x > Thre_tb, diff_indx_nb)

    if (size(overTh_index_nb)[1] == 0)

        if (size(overTh_index)[1] == 0)
            ClustT = time_intp[Indx_bi[1]]
            trusted_b = 1
        else
            ClustT = time_intp[Indx_bi[overTh_index[end]+1]]
            trusted_b = 1
        end
    else
        ClustT = time_intp[Indx_Not_bi[overTh_index_nb[1]]+1]
        trusted_b = 1
    end

    return ClustT*trusted_b
end





"""
Boipolar resident time with multi CTs

"""
function Bi_resiT_multiCT(N_clusT, time)

    Thre_tf = 60 # threshhold time to filter any noisy mulitipolar
    # unit: second

    Bi_counts = 1

    interp_linear = LinearInterpolation(time, N_clusT)
    time_intp = 0:1:last(time)
    N_clus_intp = interp_linear(time_intp)

    Indx_bi = findall(x -> x == 2, N_clus_intp)
    Bi_duraTT = size(Indx_bi,1) # unit: s

    Diff_indx_bi = diff(Indx_bi)
    Bi_counts = size(findall(t -> t > Thre_tf, Diff_indx_bi),1)+1

    return Bi_duraTT/Bi_counts

end

"""
Tell whether the bipolar state has changed its configuration

"""

function Whether_bi_changed(N_clus, Idx_CT)

    is_change = false

    Indx_bi = findall(x -> x == 2, N_clus)
    Bi_CT_idx = Idx_CT[:, Indx_bi]

    mean_bi_idx = mean(Bi_CT_idx, dims = 2)
    if mean_bi_idx != reshape(Bi_CT_idx[:,1], size(mean_bi_idx)[1],1)
        is_change = true
    end
    return is_change

end

"""
During last T mins,
test whether the state is stable bipolar state.
Criteria: State configuration not change during the 30mins
        with a tolerace of 60s leap time
"""
function Stable_bipokar_lastT(T,N_clus,Idx_CT,time)

    leap_t_thre = 60

    t_seconds = T*60
    time_during = time[end]-t_seconds
    indx_T = findall(x-> x > time_during, time)

    time_T = time[indx_T]
    N_clus_T = N_clus[indx_T]
    Idx_CT_T = Idx_CT[:,indx_T]

    indx_bipolar = findall(x->x == 2, N_clus_T)
    Diff_indx_bi = diff(indx_bipolar)

    indx_diff_indx_bi = findall(n_temp -> n_temp > 1,Diff_indx_bi)

    long_leap_time = false

    for _id_temp in 1:length(indx_diff_indx_bi)
        loca_indx_bi = indx_diff_indx_bi[_id_temp]
        leap_idx_start = indx_bipolar[loca_indx_bi]
        leap_idx_end = indx_bipolar[loca_indx_bi+1]

        leap_time = time_T[leap_idx_end] - time_T[leap_idx_start]
        if leap_time > leap_t_thre
            long_leap_time = true
        end
    end

    bi_idx_change = false

    if long_leap_time == false
        # print(size(Idx_CT_T))
        # print(size(indx_bipolar))
        # print(Idx_CT_T)
        # print(indx_bipolar)
        mean_bi_idx = mean(Idx_CT_T[:,indx_bipolar], dims = 2)
        if mean_bi_idx != reshape(Idx_CT_T[:,indx_bipolar][:,1], size(mean_bi_idx)[1],1)
            bi_idx_change = true
        end

    end

    if (bi_idx_change == false) && (long_leap_time == false)
        return true
    else
        return false
    end

end




"""
Separation Time and Initiation Time for 2 CTs
Input: distance, time steps
Output: Separation and initiation time
"""
function Bi_SepAndInitT(dis,time)
    Thre_mono, Thre_bi, Thre_t = 2, 12, 30 # unit: um, um, second
    Thre_mono_exclude = 2.5 # unit: um

    Indx_mono = findall(x -> x < Thre_mono, dis)
    Indx_bi = findall(x -> x > Thre_bi, dis)

    # Indx_last_exlucde = findall(x -> x > Thre_mono_exclude, dis[1:Indx_mono[end]])
    #
    # if Indx_last_exlucde != []
    #     Init_start = time[Indx_last_exlucde[end]]
    # else
    #     Init_start = 0
    # end

    Sep_T, Init_T = time[Indx_bi[1]], time[Indx_mono[end]]
    if Sep_T< Init_T

        # Indx_last_exlucde = findall(x -> x > Thre_mono_exclude, dis[1:Indx_mono[end]])
        #
        # if Indx_last_exlucde != []
        #     Init_start = time[Indx_last_exlucde[end]]
        # else
        #     Init_start = 0
        # end
        Indx_nextBi = findall(x -> x > Thre_bi, dis[Indx_mono[end]:end])

        if length(Indx_nextBi) > 60*20
            nextB_T = time[Indx_mono[end]:end][Indx_nextBi[1]]
            return nextB_T-Init_T, Init_T, nextB_T, Init_T
        else
            Indx_preMono = findall(x -> x < Thre_mono, dis[1:Indx_bi[1]])
            if Indx_preMono !=[]
                preMono_T = time[1:Indx_bi[1]][Indx_preMono[end]]
                return Sep_T - preMono_T,preMono_T, Sep_T, preMono_T
            end

        end
        Sep_T =  Sep_T*2 + Init_T
    end
    return Sep_T-Init_T, Init_T, Sep_T, Init_T
end
"""
Bipolar residence Time for 2 CTs
"""
function Bi_resiT(dis,time1)

    Thre_mono, Thre_t = 8, 30 # unit: um, second
    # interpolation for time step 1s
    interp_linear = LinearInterpolation(time1, dis)
    time_intp = 0:1:last(time1)
    dis_intp = interp_linear(time_intp)

    Indx_bi = findall(x -> x > Thre_mono, dis_intp)
    Bi_duraTT = size(Indx_bi,1) # unit: s

    Diff_indx_bi = diff(Indx_bi)
    Bi_counts = size(findall(t -> t > Thre_t, Diff_indx_bi),1) + 1

    return Bi_duraTT/Bi_counts
end

"""
With 2 CTs identify whether the separation is a stable state
(at least last 20mins) 1: stable bipolar; 0: unstable bipolar
"""

function BiplarStable(dis,time1)
    # Thre_mono, Thre_t = 10, 30 # unit: um, second
    # interpolation for time step 1s
    max_bi_T = 0
    max_bi_T = Bi_resiT(dis,time1)


    if max_bi_T > 30*60  # 30 mins threshold
        return (1, max_bi_T)
    else
        return (0, max_bi_T)
    end
end

## calculate the insertness for each track with multi CTs
function cal_Insertness(posi_info,CT_numb,aixs)
    insert_each = 0

    a, b, c = aixs
    T_t = length(posi_info)

    for t_temp in 1:T_t  # all the time point
        for ct_id in 1:CT_numb
            posi_t_ct_id = posi_info[t_temp][ct_id]
            # posi_pro, dis_pro = ProjDisGeo(posi_t_ct_id, axis)

            u_tar = posi_t_ct_id

            insert_temp = (u_tar[1]/a)^2 + (u_tar[2]/b)^2 + (u_tar[3]/c)^2  # outside geo >1  inside geo < 1

            insert_each = insert_each + insert_temp
        end
    end

    mean_insert = insert_each/T_t/CT_numb

    return  mean_insert
end

"""
calculate the potential energy profile in 2D geo
"""
function EnergyPotenPro2D_lxc(Para, num_CS)
    if num_CS == 2
        Potential_profile = EnergyPotentialProfile2D_lxc_2CS(Para)
    elseif num_CS == 4
        Potential_profile = EnergyPotentialProfile2D_lxc_4CS(Para)
    elseif num_CS == 0
        Potential_profile = EnergyPotentialProfile2D_lxc_0CS(Para)
    end
    return Potential_profile
end


function EnergyPotentialProfile2D_lxc_0CS(Para)
    x_1, x_2, x_3, x1_fix, k_a, k_r, l_a, l_r, p_l, k_in, k_out, a, b, c = Para
    Potential_profile = zeros(length(x_1),length(x_2))

    for n_x1 in 1:length(x_1)
        for n_x2 in 1:length(x_2)

            x_input = [x_1[n_x1], x_2[n_x2], x_3]

            x_input = [x_1[n_x1], x_2[n_x2], 0.0]
            (x_o, dis_o) = ProjDisGeo(x_input,[a,b,c])
            if  ( (x_input[1]/a)^2 + (x_input[2]/b)^2 + (x_input[3]/c)^2 )> 1
                P_ra = RadialEnergy(dis_o,k_in)
            else
                if isnan(dis_o)
                    P_ra = RadialEnergy(c,k_out)
                else
                    P_ra = RadialEnergy(dis_o,k_out)
                end

            end

            Potential_profile[n_x1,n_x2] = P_ra

        end
    end
    l_x1 = length(x_1)
    l_x2 = length(x_2)
    Potential_profile[l_x1÷2+1,l_x2÷2+1] = (Potential_profile[l_x1÷2,l_x2÷2+1]+Potential_profile[l_x1÷2+2,l_x2÷2+1])/2
    # Potential_profile = Potential_profile - minimum(Potential_profile)*ones(l_x1,l_x2)
    ## artifitial change
    mid_value = Potential_profile[l_x1÷2+1,l_x2÷2+1]
    for i_ in 1:l_x1
        for j_ in 1:l_x2
            if Potential_profile[i_, j_] > mid_value
                Potential_profile[i_, j_] = mid_value
            end
        end
    end



    return Potential_profile
end





function EnergyPotentialProfile2D_lxc_2CS(Para)
    x_1, x_2, x_3, x1_fix, k_a, k_r, l_a, l_r, p_l, k_in, k_out, a, b, c = Para
    Potential_profile = zeros(length(x_1),length(x_2))

    for n_x1 in 1:length(x_1)
        for n_x2 in 1:length(x_2)

            x_input = [x_1[n_x1], x_2[n_x2], x_3]
            dis = euclidean(x_input, x1_fix)
            P_ct = InterCTenergy(dis,k_a,k_r,l_a,l_r,p_l)

            x_input = [x_1[n_x1], x_2[n_x2], 0.0]
            (x_o, dis_o) = ProjDisGeo(x_input,[a,b,c])
            if  ( (x_input[1]/a)^2 + (x_input[2]/b)^2 + (x_input[3]/c)^2 )> 1
                P_ra = RadialEnergy(dis_o,k_in)
            else
                if isnan(dis_o)
                    P_ra = RadialEnergy(c,k_out)
                else
                    P_ra = RadialEnergy(dis_o,k_out)
                end

            end

            Potential_profile[n_x1,n_x2] = P_ct + P_ra

        end
    end
    l_x1 = length(x_1)
    l_x2 = length(x_2)
    Potential_profile[l_x1÷2+1,l_x2÷2+1] = (Potential_profile[l_x1÷2,l_x2÷2+1]+Potential_profile[l_x1÷2+2,l_x2÷2+1])/2
    # Potential_profile = Potential_profile - minimum(Potential_profile)*ones(l_x1,l_x2)

    ## artifitial change
    mid_value = Potential_profile[l_x1÷2+1,l_x2÷2+1]
    for i_ in 1:l_x1
        for j_ in 1:l_x2
            if Potential_profile[i_, j_] > mid_value
                Potential_profile[i_, j_] = mid_value
            end
        end
    end




    return Potential_profile
end


function EnergyPotentialProfile2D_lxc_4CS(Para)
    x_1, x_2, x_3, xs_fix, k_a, k_r, l_a, l_r, p_l, k_in, k_out, a, b, c = Para
    Potential_profile = zeros(length(x_1),length(x_2))

    num_CS = length(xs_fix)

    for n_x1 in 1:length(x_1)
        for n_x2 in 1:length(x_2)

            x_input = [x_1[n_x1], x_2[n_x2], x_3]
            P_ct = 0.0

            for i_ in 1:num_CS
                x1_fix = xs_fix[i_]
                dis = euclidean(x_input, x1_fix)
                P_ct = P_ct + InterCTenergy(dis,k_a,k_r,l_a,l_r,p_l)
            end


            x_input = [x_1[n_x1], x_2[n_x2], 0.0]
            (x_o, dis_o) = ProjDisGeo(x_input,[a,b,c])
            if  ( (x_input[1]/a)^2 + (x_input[2]/b)^2 + (x_input[3]/c)^2 )> 1
                P_ra = RadialEnergy(dis_o,k_in)
            else
                if isnan(dis_o)
                    P_ra = RadialEnergy(c,k_out)
                else
                    P_ra = RadialEnergy(dis_o,k_out)
                end

            end

            Potential_profile[n_x1,n_x2] = P_ct + P_ra

        end
    end
    l_x1 = length(x_1)
    l_x2 = length(x_2)
    Potential_profile[l_x1÷2+1,l_x2÷2+1] = (Potential_profile[l_x1÷2,l_x2÷2+1]+Potential_profile[l_x1÷2+2,l_x2÷2+1])/2
    # Potential_profile = Potential_profile - minimum(Potential_profile)*ones(l_x1,l_x2)
    ## artifitial change
    mid_value = Potential_profile[l_x1÷2+1,l_x2÷2+1]
    for i_ in 1:l_x1
        for j_ in 1:l_x2
            if Potential_profile[i_, j_] > mid_value
                Potential_profile[i_, j_] = mid_value
            end
        end
    end


    return Potential_profile
end
