from numba import jit

@jit(nopython=True, nogil=True)
def rmodes(v_recon_orig,v_recon_sim,max_no_segments,first_segment,inc_segments,len_time_segments,k_bottom_w,v_sim_seq,v_fit_seq,vert_mode,merid_mode):

    for iseg in range ( max_no_segments ): 
        it_seg_start = first_segment + iseg * inc_segments
        for it in range ( len_time_segments ) : 
            itime = it_seg_start + it 
            for kz in range ( k_bottom_w ) :  
                amp =  v_sim_seq[iseg,it] * vert_mode[kz] #/ ( gdef * rho_0 )  # see documentation eqn (9) for denominator
                v_recon_orig[itime, kz, :] = v_recon_orig[itime, kz, :] + amp * merid_mode[:] #* dy_tilde
                amp =  v_fit_seq[iseg,it] * vert_mode[kz] #/ ( gdef * rho_0 )  # see documentation eqn (9) for denominator
                v_recon_sim[itime, kz, :] = v_recon_sim[itime, kz, :] + amp * merid_mode[:] #* dy_tilde

    return v_recon_orig,v_recon_sim

@jit(nopython=True, nogil=True)
def calc_moc(moc,mocR,mocRt,v,v_recon_orig,v_recon_sim,e3t,k_bottom_w):

    for tt in range (v.shape[0]):
        for jy in range (v.shape[2]):
            for jk in range (v.shape[1]-1):
                moc[tt,jk+1,jy] = moc[tt,jk,jy] + v[tt,jk,jy]*e3t[jk]
                if jk < k_bottom_w-1:
                   mocR[tt,jk+1,jy]  = mocR[tt,jk,jy] + v_recon_orig[tt,jk,jy]*e3t[jk]
                   mocRt[tt,jk+1,jy] = mocRt[tt,jk,jy] + v_recon_sim[tt,jk,jy]*e3t[jk]

    return moc/1e6,mocR/1e6,mocRt/1e6
