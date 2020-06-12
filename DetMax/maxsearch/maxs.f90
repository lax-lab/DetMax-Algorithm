program main
    use omp_lib
    use algo
    use dataload
    use DFWIN
    implicit doubleprecision(a-h,o-z)
    doubleprecision time_begin,time_end1,time_end2
    integer(8) itotal,ic,multi,n_Gauge(2)
    integer(kind = 4) number_of_threads 
    double precision a_matrix(4,4),b_matrix(4,4),c_matrix(4,4),pos(3,3)
    character(len=200)::arg,ctemp,ematfile
    integer ::ISTAT,num_Gauge,PreC
    common ISTAT,num_Gauge,PreC
    call CPU_TIME(time_begin)
    narg=IARGC()
    if(narg.gt.0)then
        call getarg(1,arg)
        rptfile = arg
        call getarg(2,arg)
        nodefile = arg
        call getarg(3,arg)
        outfile = arg
        call getarg(4,arg)
        ematfile = arg
        call getarg(5,arg)
        ctemp = arg
        read(ctemp,*)ISTAT
        call getarg(6,arg)
        ctemp = arg
        read(ctemp,*)n_Gauge(1)
        call getarg(7,arg)
        ctemp = arg
        read(ctemp,*)n_Gauge(2)
        call getarg(8,arg)
        ctemp = arg
        read(ctemp,*)iTerate
        call getarg(9,arg)
        ctemp = arg
        read(ctemp,*)ipre
        call getarg(10,arg)
        ctemp = arg
        read(ctemp,*)number_of_threads
    endif
    PreC = 180/ipre
    
    write(*,"(A71)")" +--------------------------------------------------------------------+"
    write(*,"(A2,T22,A32,T71,A1)")" |","Command Parameters Input for Run","|"
    write(*,"(A71)")" +--------------------------------------------------------------------+"
    write(*,"(A71)")" +--------------------------------------------------------------------+"
    write(*,"(A2,T10,A13,T48,A1,T53,A4,T60,A1,T65,A4,T71,A1)")" |","Variable Name","|","Arg1","|","Arg2","|"
    write(*,"(A71)")" +--------------------------------------------------------------------+"
    if(n_Gauge(1).ge.1)then
        write(*,"(A2,T5,A16,T48,A1,T53,I2,T60,A1,T66,A1,T71,A1)")" |","Min Guage Number","|",n_Gauge(1),"|","-","|"
        write(*,"(A2,T5,A16,T48,A1,T53,I2,T60,A1,T66,A1,T71,A1)")" |","Max Guage Number","|",n_Gauge(2),"|","-","|"
    else
        write(*,"(A2,T5,A12,T48,A1,T53,I2,T60,A1,T66,A1,T71,A1)")" |","Guage Number","|",n_Gauge(1),"|","-","|"
    endif
    write(*,"(A2,T5,A23,T48,A1,T53,I2,T60,A1,T66,A1,T71,A1)")" |","DetMax Iteration Number","|",iTerate,"|","-","|"
    write(*,"(A2,T5,A11,T48,A1,T53,I2,T60,A1,T66,A1,T71,A1)")" |","Load Number","|",ISTAT,"|","-","|"
    write(*,"(A2,T5,A21,T48,A1,T53,I2,T60,A1,T66,A1,T71,A1)")" |","Angle Search Accuracy","|",ipre,"|","-","|"
    write(*,"(A2,T5,A19,T48,A1,T53,I2,T60,A1,T66,A1,T71,A1)")" |","Number of Processor","|",number_of_threads,"|","-","|"
    write(*,"(A2,T5,A21T48,A1,T52,f5.3,T60,A1,T66,A1,T71,A1)")" |","Convergence Criterion","|",0.001,"|","-","|"
    write(*,"(A71)")" +--------------------------------------------------------------------+"
    
    call dataframe() 
    call search_matrix()
    open(unit = 10,file = outfile)       
    open(unit = 11,file = ematfile)    
    call OMP_set_num_threads(number_of_threads)
    if(n_Gauge(2).lt.1)then
        num_Gauge = n_Gauge(1) 
        allocate(findmax(number_of_threads),detsel(number_of_threads))
        allocate(getmax%icho(num_Gauge),getmax%cita(num_Gauge))
        allocate(getmax%F_matrix(num_Gauge,ISTAT),getmax%sM_inv(ISTAT,ISTAT))
        do iall=1,number_of_threads
            allocate(findmax(iall)%icho(num_Gauge),findmax(iall)%cita(num_Gauge))
            allocate(findmax(iall)%F_matrix(num_Gauge,ISTAT),findmax(iall)%sM_inv(ISTAT,ISTAT))
            allocate(detsel(iall)%icho(num_Gauge),detsel(iall)%cita(num_Gauge))
            allocate(detsel(iall)%F_matrix(num_Gauge,ISTAT),detsel(iall)%sM_inv(ISTAT,ISTAT))
        enddo
        getmax%detM = 0.0
        
        !$OMP PARALLEL DEFAULT(PRIVATE),SHARED(number_of_threads,findmax, detsel,search_list) 
        !$OMP DO        
        do i=1,number_of_threads            
            iProc = OMP_get_thread_num()+1
            call detmax(iProc)
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        do i=1,number_of_threads
            if(getmax%detM.lt.findmax(i)%detM)then
                getmax = findmax(i)
            endif
        enddo
        do i=1,num_Gauge
            write(11,"(I3,2X,I7,2X,I5,2X,20(ES11.4,2X))")i,search_list(getmax%icho(i))%id,getmax%cita(i),(getmax%F_matrix(i,j),j=1,ISTAT)
        enddo
        write(10,*)"Detmax:",getmax%detM
        write(10,"(10(I7,5X))")(search_list(getmax%icho(i))%id,i=1,num_Gauge)
        write(10,"(10(I7,5X))")(getmax%cita(i),i=1,num_Gauge)
        write(10,*)"|-----------------------Cordinate-------------------------|"        
        do i=1,num_Gauge
            write(10,*)"Guage",i,"£º"
            no_ele = search_list(getmax%icho(i))%id
            theta = getmax%cita(i)
            call Roderick_Trans(no_ele,theta,pos)
            do k=1,3
                write(11,*)(pos(k,j),j=1,3)
            enddo
        enddo
        deallocate(getmax%icho,getmax%cita,getmax%F_matrix,getmax%sM_inv)
        do iall=1,number_of_threads
            deallocate(findmax(iall)%icho,detsel(iall)%icho)
            deallocate(findmax(iall)%cita,detsel(iall)%cita)
            deallocate(findmax(iall)%F_matrix,detsel(iall)%F_matrix)
            deallocate(findmax(iall)%sM_inv,detsel(iall)%sM_inv)
        enddo
        deallocate(findmax,detsel)
        
    else
        do ig =n_Gauge(1),n_Gauge(2)
            num_Gauge = ig
            allocate(findmax(number_of_threads),detsel(number_of_threads))
            allocate(getmax%icho(num_Gauge),getmax%cita(num_Gauge))
            allocate(getmax%F_matrix(num_Gauge,ISTAT),getmax%sM_inv(ISTAT,ISTAT))
            do iall=1,number_of_threads
                allocate(findmax(iall)%icho(num_Gauge),findmax(iall)%cita(num_Gauge))
                allocate(findmax(iall)%F_matrix(num_Gauge,ISTAT),findmax(iall)%sM_inv(ISTAT,ISTAT))
                allocate(detsel(iall)%icho(num_Gauge),detsel(iall)%cita(num_Gauge))
                allocate(detsel(iall)%F_matrix(num_Gauge,ISTAT),detsel(iall)%sM_inv(ISTAT,ISTAT))
            enddo
            getmax%detM = 0.0
            !$OMP PARALLEL DEFAULT(PRIVATE),SHARED(number_of_threads,findmax, detsel,search_list) 
            !$OMP DO        
            do i=1,number_of_threads
                iProc = OMP_get_thread_num()+1
                call detmax(iProc)
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            do i=1,number_of_threads
                if(getmax%detM.lt.findmax(i)%detM)then
                    getmax = findmax(i)
                endif
            enddo 
            do i=1,num_Gauge
                write(11,"(I3,2X,I7,2X,I5,2X,20(ES11.4,2X))")i,search_list(getmax%icho(i))%id,getmax%cita(i),(getmax%F_matrix(i,j),j=1,ISTAT)
            enddo
            write(10,*)"Detmax:",getmax%detM
            write(10,"(10(I7,5X))")(search_list(getmax%icho(i))%id,i=1,num_Gauge)
            write(10,"(10(I7,5X))")(getmax%cita(i),i=1,num_Gauge)
            write(10,*)"|-----------------------Cordinate-------------------------|"               
            do i=1,num_Gauge
                write(10,*)"Guage",i,"£º"
                no_ele = search_list(getmax%icho(i))%id
                theta = getmax%cita(i)
                call Roderick_Trans(no_ele,theta,pos)
                do k=1,3
                    write(11,*)(pos(k,j),j=1,3)
                enddo
            enddo
            deallocate(getmax%icho,getmax%cita,getmax%F_matrix,getmax%sM_inv)
            do iall=1,number_of_threads
                deallocate(findmax(iall)%icho,detsel(iall)%icho)
                deallocate(findmax(iall)%cita,detsel(iall)%cita)
                deallocate(findmax(iall)%F_matrix,detsel(iall)%F_matrix)
                deallocate(findmax(iall)%sM_inv,detsel(iall)%sM_inv)
            enddo
            deallocate(findmax,detsel)          
        enddo
    endif
    close(11)
    call CPU_TIME(time_end1)
   
    end program

    
    
