import numpy as np
from scipy.special import comb,perm

core_orbit = 4
valence_orbit = 4
g = 0.5


def v(p,q,r,s):
    if ( (int(p/2) == int(q/2)) and (int(r/2) == int(s/2)) and ( p !=q ) and ( r!= s) ):
        if ( (p < q) and (r < s)):
            return -g/2
        elif ( (p < q) and (r > s)):
            return g/2
        elif ( (p > q) and (r < s)):
            return g/2
        elif ( (p > q) and (r > s)):
            return -g/2
        else:
            return 0
    else: 
        return 0 


def spe(p,q):
    d = 1
    if (p == q):
        return int(p/2)*d
    else:
        return 0

def f_pq(p,q):
    sum_v_piqi = 0
    for i in range (0,core_orbit_num):
        sum_v_piqi = sum_v_piqi + v(p,i,q,i)
    f_pq = single_particle(p,q)+sum_v_piqi
    return f_pq 



def CCD_main(orbit_ij,orbit_ab):
    tot_orbit = orbit_ij + orbit_ab
    H_abij = np.zeros((orbit_ab,orbit_ab,orbit_ij,orbit_ij))
    t_abij = np.zeros((orbit_ab,orbit_ab,orbit_ij,orbit_ij))
    f_pq_matrix = np.zeros((tot_orbit,tot_orbit))
    v_pqrs = np.zeros((tot_orbit,tot_orbit,tot_orbit,tot_orbit))
    for p in range(0,tot_orbit):
        for q in range(0,tot_orbit):
            for r in range(0,tot_orbit):
                for s in range(0,tot_orbit):
                    v_pqrs[p,q,r,s] = v(p,q,r,s)

    #for a in range(0,orbit_ab):
    #    for b in range(0,orbit_ab):
    #        for i in range(0,orbit_ij):
    #            for j in range(0,orbit_ij):
    #                t_abij [a,b,i,j] = v(a,b,i,j)/(spe(a,a)+spe(b,b)-spe(i,i)-spe(j,j)) 
#                h_1 = v(a,b,i,j)
    #                h_2 = (single_particle(a,a)+single_particle(b,b)+single_particle(c,c)+single_particle(d,d))*t_abij[a,b,i,j]
    #                h_3 = 0.5*
    #
      
    #h1 = 0.5*np.einsum('abcd,cdij->abij',v[ orbit_ij:tot_orbit, orbit_ij:tot_orbit,orbit_ij:tot_orbit,orbit_ij:tot_orbit ],t_abij)
    h1 = v_pqrs[orbit_ij:tot_orbit,orbit_ij:tot_orbit,0:orbit_ij,0:orbit_ij]
    print ("h1="+str(h1))
    return h1[0,0,0,0]


h1 = CCD_main(core_orbit,valence_orbit)
print("h1="+str(h1))
g
#for p in range (0,4):
#    for q in range (0,4):
#        print ("f"+str(p)+str(q)+"="+str(f_pq(p,q)))

#for p in range(0,4):
#    for q in range(0,4):
#        for r in range(0,4):
#            for s in range(0,4):
#                print ("v "+str(p)+str(q)+str(r)+str(s)+"="+str(v(p,q,r,s)))
