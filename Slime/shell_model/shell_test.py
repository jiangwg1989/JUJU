import numpy as np
from scipy.special import comb,perm


a = 0b11000

#print(a)
a = 10
b = 13
print("np="+str(np.bitwise_and(a,b)))


#
#   Set up all the possible configurations
#

def set_config(tot_orbit,particle_num):
    tot_config_num = int(comb(tot_orbit,particle_num))
    config_all = np.zeros(tot_config_num)
    
    config_start   =  0
    for i in range(0,particle_num):
        config_start = config_start + 2**i

    config_end = np.left_shift( int(config_start),tot_orbit-particle_num)
    print("config_start="+str(bin(config_start)))
    print("config_end="+str(bin(config_end)))

    config_all[0] = int(config_start)
    config_temp = int(config_start)
     
    #print("1count="+str(bin(config_temp+2).count("1")))
    #input()
    i = 1
    while (1):
        config_temp = config_temp + 1
        if ( bin(config_temp).count("1") == int(particle_num) ):
            config_all [i] = config_temp 
            i = i + 1
        if (i == tot_config_num):
            break
#    for i in range(0,len(config_all)):
#        print(bin(int(config_all[i])))
#        u  = int(config_all[i-1])
#        print("u="+str(bin(u)))
#        print("-u="+str(bin(-u)))
#        #u   = 15
#        ur = np.bitwise_and(u,-u)
#        print("ur="+str(bin(ur)))
#        nh = u + ur
#        print("uh="+str(bin(nh)))
#        n0 = np.bitwise_xor(u,nh)
#        print("u0="+str(bin(n0)))
#        n0 = int(n0/ur)
#        print ("n0="+str(bin(n0)))
#        #n0 = np.left_shift(n0,2)
#        print ("n0="+str(bin(n0)))
#        u  = np.bitwise_or(nh,n0)
#        config_all[i] = u
#        print("lalala="+str(bin(u)))
    #print(bin(config_all))
    #config_end     =
    return tot_config_num,config_all    

def pairing_config(config_all):
    loop1 = 0
    d = 1
    config_pairing_num = int(comb(tot_orbit/2,particle_num/2)) 
    config_pairing = np.zeros(config_pairing_num,dtype = "int")
    spe_sum = np.zeros(config_pairing_num)
    for loop2 in range(0,len(config_all)):
        vec = config_all[loop2]
        p_flag = 1
        #print("vec="+str(bin(int(vec))))
        vec_temp = vec
        spe_sum_temp = 0
        for loop3 in range(0,particle_num):
            if (loop3 % 2 == 0):
                i = where_one(int(vec_temp))
                vec_temp = vec_temp - (1 << i )
                spe_sum_temp = spe_sum_temp + int(i/2)
                if ( i % 2 == 0 ):
                    p_flag = p_flag*1
                else:
                    p_flag = p_flag*0
                    break
                #print (i)
            else:
                j = where_one(int(vec_temp))
                vec_temp = vec_temp - (1 << j )
                spe_sum_temp = spe_sum_temp + int(j/2)
                if ( j - i == 1 ):
                    p_flag = p_flag * 1
                else:
                    p_flag = p_flag * 0
                    break
                #print (j)
        if (p_flag == 1):
            config_pairing[loop1] = vec
            spe_sum [loop1] = spe_sum_temp
            print ("config_pair="+str(bin(int(vec))))
            print ("spe_sum_temp ="+str(spe_sum_temp))
            loop1 = loop1 + 1
    if (loop1 != config_pairing_num):
        print ("wtf: config_pairing_num error")
    print ("config_pairing = "+str(config_pairing)) 
    print ("spe_sum ="+str(spe_sum))
    return spe_sum,config_pairing_num,config_pairing

def matrix_element(a,b,i,j):
    # paring model 
    g = -1
    if ( (a % 2 == 0) and ( b-a == 1) and (i%2 ==0) and (j-i == 1)):
        return g/2
    else:
        return 0 


def where_one(int_num):
    bin_str = bin(int_num).replace("0b",'')
    for i in range(256):
        if (1<<i) & int_num:
            return i


tot_orbit = 8
particle_num = 4


config_all = np.zeros(int(comb(tot_orbit,particle_num)),dtype = "int")
config_pairing = np.zeros(int(comb(tot_orbit/2,particle_num/2)),dtype = "int")

config_pairing_num = int(comb(tot_orbit/2,particle_num/2)) 
spe_sum = np.zeros(config_pairing_num)
tot_config_num,config_all = set_config(tot_orbit,particle_num)
spe_sum,pair_config_num,config_pairing = pairing_config(config_all)

print ("tot_config_num="+str(tot_config_num))
print ("pair_config_num="+str(pair_config_num))
print ("spe_sum:"+str(spe_sum))
#input()



#
# calculate the H matrix elements
#
H_matrix_elements =  np.zeros((pair_config_num,pair_config_num))

for loop1 in range(0,pair_config_num):
    for loop2 in range(loop1,pair_config_num):
        #print ("loop1 2="+str(bin(int(config_all[loop1])))+","+str(bin(int(config_all[loop2]))))
        right_vec = int(config_pairing[loop1])
        left_vec  = int(config_pairing[loop2])
        dif_bit_num = bin(np.bitwise_xor(right_vec,left_vec)).count("1")/2        
#        
        if (dif_bit_num == 2):
            right_left_and = np.bitwise_and(right_vec,left_vec)            
            bit_ij = np.bitwise_xor(right_left_and,right_vec)
            bit_ab = np.bitwise_xor(right_left_and,left_vec)
            #print ("right_vec="+str(bin(right_vec)))
            #print ("left_vec ="+str(bin(left_vec)))
            #print ("bit_ab="+str(bin(bit_ab)))
            #print ("bit_ij="+str(bin(bit_ij)))
            i = where_one(bit_ij)
            bit_i = 1 << i 
            bit_j = bit_ij - bit_i
            j = where_one(bit_j)
            one_between_ij = -1 
            for loop3 in range(i,j):
                if (1<<loop3) & right_vec:
                    one_between_ij = one_between_ij + 1
            #print ("num_one_ij = "+str(one_between_ij))

            #print ("i="+str(i)+" j="+str(j))
            if (i >= j):
                print("wtf: i >= j ")
            
            a = where_one(bit_ab)
            bit_a = 1 << a 
            bit_b = bit_ab - bit_a
            b = where_one(bit_b)
            one_between_ab = -1 
            for loop3 in range(a,b):
                if (1<<loop3) & left_vec:
                    one_between_ab = one_between_ab + 1
            #print ("num_one_ab = "+str(one_between_ab))

            #print ("a="+str(a)+"  b="+str(b))
            if (a >= b):
                print("wtf: a >= b ")
            if ( (one_between_ab + one_between_ij ) % 2 == 0 ):
                phase = 1
            else:
                phase = -1
            #print ("one_between="+str(one_between_ab+one_between_ij)+" phase ="+str(phase))
            H_matrix_elements[loop1,loop2] = phase * matrix_element(a,b,i,j)
            
        if (dif_bit_num == 1):
            right_left_and = np.bitwise_and(right_vec,left_vec)            
            bit_i = np.bitwise_xor(right_left_and,right_vec)
            bit_a = np.bitwise_xor(right_left_and,left_vec)
            bit_rest = right_vec - bit_i

            #print ("right_vec="+str(bin(right_vec)))
            #print ("left_vec="+str(bin(left_vec)))
            #a = where_one()
            #b = where_one()
            i = where_one(bit_i)
            a = where_one(bit_a)
            bit_rest_temp = bit_rest
            for loop3 in range(0, particle_num - 1):
                j = where_one(bit_rest_temp)
                b = j
                bit_j = 1 << j
                bit_b = bit_j
                bit_rest_temp = bit_rest_temp - bit_j                
                #print ("bit_a="+str(bin(bit_a)))
                #print ("bit_b="+str(bin(bit_b)))
                #print ("bit_i="+str(bin(bit_i)))
                #print ("bit_j="+str(bin(bit_j)))
 
                aa = a
                bb = b
                ii = i
                jj = j
                # two "a+" or two "a" equals zero
                if (  (a ==b) or (i ==j) ):  
                    continue
                # make sure a < b and i < j
                if ( a > b ):
                    aa,bb = b, a
                if ( i > j ):
                    ii,jj = j, i

                one_between_aabb = 0 
                for loop4 in range(aa,bb-1):
                    if (1<<loop4 + 1) & bit_rest:
                        one_between_aabb = one_between_aabb + 1
                #print ("num_one_ab = "+str(one_between_aabb))
                one_between_iijj = 0 
                for loop4 in range(ii,jj-1):
                    if (1<<loop4 + 1 ) & bit_rest:
                        one_between_iijj = one_between_iijj + 1
                #print ("num_one_ij = "+str(one_between_iijj))
                if ( (one_between_aabb + one_between_iijj ) % 2 == 0 ):
                    phase = 1
                else:
                    phase = -1
                H_matrix_elements[loop1,loop2]= H_matrix_elements[loop1,loop2] + phase * matrix_element(aa,bb,ii,jj)
               
                #print ("ab ij = "+str(aa)+" "+str(bb)+" "+str(ii)+" "+str(jj))  
                #print ("one_between="+str(one_between_aabb+one_between_iijj)+" phase ="+str(phase))
            
             
        if (dif_bit_num == 0):
            if ( right_vec != left_vec ):
                print("wtf: right!= left")
            occupy_orbit = np.zeros(particle_num,dtype = "int")
            vec_temp = right_vec
            for loop3 in range(0,particle_num):
                occupy_orbit[loop3] = where_one(vec_temp)
                vec_temp = vec_temp - ( 1 << occupy_orbit[loop3] )
            #print (bin(right_vec))
            #print (occupy_orbit)        
            for loop3 in range(0,particle_num):
                for loop4 in range(loop3+1,particle_num):
                    a = occupy_orbit[loop3]
                    b = occupy_orbit[loop4]
                    i = a 
                    j = b
                    #print ("ab ij = "+str(a)+" "+str(b)+" "+str(i)+" "+str(j))  
                    H_matrix_elements[loop1,loop2]= H_matrix_elements[loop1,loop2] + matrix_element(a,b,i,j)
            #print("right_vec="+str(bin(right_vec)))
            #print("H_elements="+str(H_matrix_elements[loop1,loop2]))      
for loop1 in range(0,pair_config_num):
    for loop2 in range(loop1+1,pair_config_num):
        H_matrix_elements[loop2,loop1] = H_matrix_elements[loop1,loop2]

#
# including single particle energy
#
for loop1 in range(0,pair_config_num):
    H_matrix_elements[loop1,loop1] = H_matrix_elements[loop1,loop1] + spe_sum[loop1] 

print (H_matrix_elements)
             
#for a in range(0,tot_orbit):             
#    for b in range(a,tot_orbit):             
#        for i in range(0,tot_orbit):             
#            for j in range(i,tot_orbit):
#                test = matrix_element(a,b,i,j)
#                if (test != 0):
#                    print ("ab ij = "+str(a)+" "+str(b)+" "+str(i)+" "+str(j))
#                    print(test)

#
# diagonalization
#
eigen_value,eigen_vec = np.linalg.eig(H_matrix_elements)
gs_energy = np.min(eigen_value)
print("gs_energy = ", gs_energy)

