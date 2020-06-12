	module algo
    use dataload
    contains
        !求矩阵相乘
        subroutine brmul(a_matrix,b_matrix,intm,intn,intk,c_matrix)
        implicit double precision(a-h,o-z)
        double precision a_matrix(intm,intn),b_matrix(intn,intk)
        double precision c_matrix(intm,intk)
        do i=1,intm
            do j=1,intk
                c_matrix(i,j)=0.0
                do k=1,intn
                    c_matrix(i,j)=c_matrix(i,j)+a_matrix(i,k)*b_matrix(k,j)
                enddo
            enddo
        enddo
        end subroutine brmul
!        求取行列式的值
!       a_matrix的值被修改，需要定义
        subroutine bsdet(t_matrix,intn,det)
        implicit double precision(a-h,o-z)
        double precision a_matrix(intn,intn),t_matrix(intn,intn)
        double precision det,fvalue,dvalue,qvalue
        a_matrix = t_matrix
        fvalue =1.0
        det =1.0
        do k=1,intn-1
            qvalue = 0.0
            do i=k,intn
                do j=k,intn
                    if(abs(a_matrix(i,j)).gt.qvalue)then
                        qvalue = abs(a_matrix(i,j))
                        is=i
                        js=j
                    endif
                enddo
            enddo
            if(qvalue.eq.0.0)then
                det=0.0
                return
            endif
            if(is.ne.k)then
                fvalue = -fvalue
                do j=k,intn
                    dvalue = a_matrix(k,j)
                    a_matrix(k,j) = a_matrix(is,j)
                    a_matrix(is,j) = dvalue
                enddo
            endif
            if(js.ne.k)then
                fvalue = -fvalue
                do i=k,intn
                    dvalue = a_matrix(i,js)
                    a_matrix(i,js) = a_matrix(i,k)
                    a_matrix(i,k) = dvalue
                enddo
            endif
            det = det*a_matrix(k,k)
            do i=k+1,intn
                dvalue = a_matrix(i,k)/a_matrix(k,k)
                do j=k+1,intn
                    a_matrix(i,j) = a_matrix(i,j)-dvalue*a_matrix(k,j)
                enddo
            enddo
        enddo
        det = fvalue*det*a_matrix(intn,intn)
        end subroutine bsdet  
        
        subroutine detsolve(a_matrix,intn,det)
        implicit double precision(a-h,o-z)
        double precision a_matrix(intn,intn),l_matrix(intn,intn),u_matrix(intn,intn)        
        call crout(a_matrix,l_matrix,u_matrix,intn)
        det = 1.0
        do i=1,intn
            det = det*l_matrix(i,i)
        enddo
        end subroutine detsolve
        !Lu crout 分解
        subroutine crout(a_matrix,l_matrix,u_matrix,intn)
        double precision a_matrix(intn,intn),l_matrix(intn,intn),u_matrix(intn,intn)
        l_matrix(:,1) = a_matrix(:,1)
        u_matrix(1,:)=a_matrix(1,:)/l_matrix(1,1)
        do k=2,intn
            do i=k,intn
                s=0.0
                do m =1,k-1
                    s=s+l_matrix(i,m)*u_matrix(m,k)
                enddo
                l_matrix(i,k) = a_matrix(i,k)-s
            enddo
            do j=k+1,intn
                s=0.0
                do m=1,k-1
                    s=s+l_matrix(k,m)*u_matrix(m,j)
                enddo
                u_matrix(k,j) = (a_matrix(k,j)-s)/l_matrix(k,k)
            enddo
            u_matrix(k,k)=1.0
        enddo
        end subroutine crout        
        !求取实对称矩阵的逆
        !调用说明 输入的matrix_A会被修改，因此在调用前需采用变量赋值方式传入
        !temp = matrix_A
        !call bssgj(temp,2)
        subroutine bssgj(matrix_A,intn)
        implicit double precision(a-h,o-z)
        double precision matrix_A(intn,intn),vector_B(intn)
        double precision value_w,value_g
        intl =1
        do k=1,intn
            m=intn-k+1
            value_w=matrix_A(1,1)
            !注意防止出现程序将无穷小值认为i是零
            if(value_w.eq.0.0)then
                intl=0
                !write(*,*)">>>Warning,Matrix inverse Failed!"
                !write(*,*)">>>Renitializing..!"
                return
            endif
            do i=2,intn
                value_g=matrix_A(i,1)
                vector_B(i)=value_g/value_w
                if(i.le.m)vector_B(i)=-vector_B(i)
                do j=2,i
                    matrix_A(i-1,j-1)=matrix_A(i,j)+value_g*vector_B(j)
                enddo
            enddo
            matrix_A(intn,intn)=1.0/value_w
            do i=2,intn
                matrix_A(intn,i-1)=vector_B(i)
            enddo
        enddo
        do i=1,intn-1
            do j=i+1,intn
                matrix_A(i,j)=matrix_A(j,i)
            enddo
        enddo        
        end subroutine bssgj
        !向量叉乘
        subroutine cross(vector_A,vector_B,vector_C)
        implicit double precision(a-h,o-z)
        double precision vector_A(3),vector_B(3)
        double precision vector_C(3)
        vector_C(1) = vector_A(2)*vector_B(3)-vector_A(3)*vector_B(2)
        vector_C(2) = vector_A(3)*vector_B(1)-vector_A(1)*vector_B(3)
        vector_C(3) = vector_A(1)*vector_B(2)-vector_A(2)*vector_B(1)
        end subroutine cross
        !D-optimal得到对应的贴片位置和贴片参数后（局部坐标系），通过将全局坐标系X轴或Z轴投影到ele平面
        !作为ele的x方向，右手系90°为y方向，旋转角度方向通过罗德里格变换获得
        !空间向量绕ele平面法向旋转对应角度
        !如果在定义截面属性时指定了壳单元的局部坐标系，则应力、应变和截面力分量的方向都是基于此局部坐标系。
        !如果没有指定单元的局部坐标系，则采用默认的局部坐标系方向，其确定方法如下：
        !将全局坐标系的x方向在壳单元面上的投影作为局部坐标系的1方向。如果全局坐标系的x方向与壳单元面法线方向的夹角≤0.1°，
        !则将全局坐标系的z方向在壳单元面上的投影作为局部坐标第的1方向（如下图右）。
        !input：
        !单元号，D-optimal角度
        subroutine Roderick_Trans(no_ele,theta,pos)
        implicit double precision(a-h,o-z)
        integer ipos_node(3)
        double precision matrix_a(3,2)
        double precision vector_A(3),vector_B(3),vector_C(3)
        double precision cord(3),projx(3),C_norm(3),axis_y(3)
        double precision ata(2,2),step1(3,2),step2(3,3)
        double precision res_vec(3),pos(3,3)
        ifind = -1
        itotal = size(ele,dim=2)
        iter = 1
        !找到单元号所在的位置
        do while(ifind.lt.0)
            minus = abs(ele(1,iter)%id-no_ele)
            if(minus.gt.0)then
                iter=iter+1
            else
                ifind = 1
                ipos_ele = iter
            endif
        enddo
        !找到单元号对应node的位置
        itotal=size(nod)
        do i=1,3
            ifind = -1
            no_nod = ele(1,ipos_ele)%node(i)
            iter = 1
            do while(ifind.lt.0)
                minus = abs(nod(iter)%id-no_nod)
                if(minus.gt.0)then
                    iter=iter+1
                else
                    ifind = 1
                    ipos_node(i)=iter
                endif
            enddo
        enddo
        !赋值基向量
        do i=1,3
            pos(i,1) = nod(ipos_node(1))%pos(i)
            vector_A(i) = nod(ipos_node(2))%pos(i)-nod(ipos_node(1))%pos(i)
            vector_B(i) = nod(ipos_node(3))%pos(i)-nod(ipos_node(2))%pos(i)
        enddo
        !赋值基矩阵        
        do i=1,3
            matrix_a(i,1)=vector_A(i)
            matrix_a(i,2)=vector_B(i)
        enddo
        ata = matmul(transpose(matrix_a),matrix_a)
        call bssgj(ata,2)
        step1 = matmul(matrix_a,ata)
        step2 = matmul(step1,transpose(matrix_a))
        !判断全局坐标系x轴与法向量夹角
        call  cross(vector_A,vector_B,vector_C)
        cord(1)=1;cord(2:3)=0
        va_norm = sqrt(vector_C(1)**2+vector_C(2)**2+vector_C(3)**2)
        do i=1,3
            C_norm(i) = vector_C(i)/va_norm
        enddo
        cri_ang = acos(dot_product(cord,C_norm))*180.0/PI
        if(cri_ang.gt.179)then
            delta = abs(cri_ang-180.0)
        else
            delta = cri_ang
        endif
        if(delta<0.1)then
            cord(1:2)=0;cord(3)=1
        endif
        !开始投影
        projx = matmul(step2,cord)
        call cross(C_norm,projx,axis_y)
        do i=1,3
            res_vec(i) = projx(i)*cos(theta*PI/180.0)+axis_y(i)*sin(theta*PI/180.0)
        enddo
        va_norm= sqrt(res_vec(1)**2+res_vec(2)**2+res_vec(3)**2)
        do i=1,3
            res_vec(i) = res_vec(i)/va_norm
        enddo
        !计算贴片方向第二点/第三点
        delta_L = 2.0
        do i=1,3
            pos(i,2) = res_vec(i)*delta_L+pos(i,1)
            pos(i,3) = res_vec(i)*delta_L+pos(i,2)
        enddo        
        end subroutine Roderick_Trans
        !建立搜索矩阵
        subroutine search_matrix()
        implicit double precision(a-h,o-z)
        integer itot
        double precision delt_cita,tr_rad
        double precision step1,step2,step3
        integer ::ISTAT,num_Gauge,PreC
        common ISTAT,num_Gauge,PreC
        !The selected groups should have some positive and some negative strain 
        !components of high magnitude under any given unit load.Each group conforms 
        !to the limitations already mentioned. Experience and engineering judgment 
        !will restrict the choice to a manageable number
        !目前试行全局范围内寻优，如果后期计算效率偏低，
        !考虑以stat1的主应变方向优化代码
        itot = size(ele,dim=2)
        delt_cita = 180.0/PreC
        tr_rad = PI/180.0
        do i = 1,itot
            search_list(i)%id = ele(1,i)%id
            search_list(i)%ichosen = 0
            do j = 1,PreC                
                search_list(i)%cita(j) = delt_cita*j
                t_cita = delt_cita*j*tr_rad
                !generate the row vector                
                do k = 1,ISTAT
                    step1 = 0.5*(ele(k,i)%strain(1)+ele(k,i)%strain(2))
                    step2 = 0.5*(ele(k,i)%strain(1)-ele(k,i)%strain(2))*cos(2.0*t_cita)
                    step3 = 0.5*ele(k,i)%strain(3)*sin(2.0*t_cita)
                    search_list(i)%eps(j,k) = step1+step2+step3
                enddo
            enddo 
        enddo
        end subroutine search_matrix 
        !随机初始化
        !方向取stat1对应的主应力方向
        !对数据list进行标定（ichosen）
        subroutine rand_init(iProc)
        implicit double precision(a-h,o-z)
        integer itot,idif,m_sita(2)
        integer,allocatable:: icho(:)
        double precision step1,step2,step3,det,comp(2)
        double precision,allocatable:: se_vec(:,:)
        double precision,allocatable:: sM_ori(:,:)
        integer ::ISTAT,num_Gauge,PreC
        common ISTAT,num_Gauge,PreC
        tr_rad = PI/180.0
        tr_sita = 180.0/PI
        search_list(:)%ichosen = 0
        itot = size(ele,dim=2)
        allocate(icho(num_Gauge))
        allocate(se_vec(2,ISTAT))
        allocate(sM_ori(ISTAT,ISTAT))
        call random_seed () 
        do i=1,num_Gauge            
            call random_number (xrand)
            num = floor(xrand*itot)+1
            if(num.gt.itot)num=num-1
            if(i.lt.2)then
                icho(i) = num
            else
                isig = 0
                do while(isig.lt.1)
                    do j=1,i-1
                        idif = abs(num - icho(j))
                        if(idif.lt.1)isig = 1
                    enddo
                    if(isig.lt.1)then
                        icho(i) = num
                        isig = 1
                    else
                        call random_number (xrand)
                        num = floor(xrand*itot)+1
                        if(num.gt.itot)num=num-1
                        isig = 0
                    endif
                enddo
            endif
            search_list(icho(i))%ichosen = 1
        enddo
        do i=1,num_Gauge
            ip = icho(i)
             tan_in = ele(1,ip)%strain(3)/(ele(1,ip)%strain(1)-ele(1,ip)%strain(2))
             msita = 0.5*tr_sita*atan(tan_in)
             if(msita.lt.0)then
                 m_sita(1) = msita+90
                 m_sita(2) = msita+180
             else
                 m_sita(1) = msita
                 m_sita(2) = msita+90
             endif
             comp(:)=0.0
             do j=1,2 
                 t_cita = m_sita(j)*tr_rad
                 do k = 1,ISTAT
                    step1 = 0.5*(ele(k,ip)%strain(1)+ele(k,ip)%strain(2))
                    step2 = 0.5*(ele(k,ip)%strain(1)-ele(k,ip)%strain(2))*cos(2.0*t_cita)
                    step3 = 0.5*ele(k,ip)%strain(3)*sin(2.0*t_cita)
                    se_vec(j,k) = step1+step2+step3
                    comp(j)=comp(j)+abs(se_vec(j,k))
                 enddo
             enddo
             if(comp(1).gt.comp(2))then
                 ic = 1
             else
                 ic = 2
             endif 
             detsel(iProc)%icho(i) =  icho(i)
             detsel(iProc)%cita(i) = m_sita(ic)
             detsel(iProc)%F_matrix(i,:) = se_vec(ic,:)
        enddo        
        !    detsel(iProc)%icho = (/1380,14876,25606,22472/)
        !    detsel(iProc)%cita = (/180,85,90,84/)
        !do i=1,num_Gauge
        !    detsel(iProc)%F_matrix(i,:)=search_list(detsel(iProc)%icho(i))%eps(detsel(iProc)%cita(i),:)
        !enddo
        sM_ori = matmul(transpose(detsel(iProc)%F_matrix),detsel(iProc)%F_matrix)        
        call bsdet(sM_ori,ISTAT,det)        
        !生成矩阵病态严重,regenerate
         if(det.lt.singular_cri)then
            det = 0.0            
        else       
            !get inverse of sM_ori, sM_ori = sM_inv
            call bssgj(sM_ori,ISTAT)
            detsel(iProc)%sM_inv = sM_ori
            iter=0
            !$OMP CRITICAL
            write(*,"(A71)")" ======================================================================"
            write(*,"(A17,T22,A3,I7,A30,I7)")"PROCESSOR","=",iProc,"OUTER LOOP ITERATION = ",iter
            write(*,"(A71)")" ----------------------------------------------------------------------"
            write(*,"(A19,T22,A3,ES11.4)")"DETERMINANT","=",det
            write(*,"(A71)")" ----------------------------------------------------------------------"
            write(*,"(A2,T12,A12,T34,A1,T50,A5,T71,A1)")"|","Guage Number","|","Angle","|"
            write(*,"(A71)")" ----------------------------------------------------------------------"
            do i=1,num_Gauge
                write(*,"(A2,T13,I7,T34,A1,T47,I7,T71,A1)")" |",search_list(detsel(iProc)%icho(i))%id,"|",detsel(iProc)%cita(i),"|"
            enddo
            write(*,"(A71)")" +-------------------------------+------------------------------------+//"
            !$OMP END CRITICAL
        endif
        detsel(iProc)%detM = det
        end subroutine rand_init 
        
        !对trapped matrix 进行行扰动
        subroutine disturb(iProc)
        implicit double precision(a-h,o-z)
        integer itot,itur(3)
        double precision,allocatable:: sM_ori(:,:) 
        integer ::ISTAT,num_Gauge,PreC
        common ISTAT,num_Gauge,PreC
        call random_seed () 
        call random_number (xrand) 
        if(xrand.gt.0.5)then
            call rand_init(iProc)
        else
            allocate(sM_ori(ISTAT,ISTAT))
            itot = size(ele,dim=2)
            call random_number (xrand) 
            isig = 0
            do while(isig.ne.1)
                itur(1) = floor(xrand*itot)+1
                if(itur(1).gt.itot)itur(1) =itur(1) -1
                if(search_list(itur(1))%ichosen.ne.1)isig = 1
                call random_number (xrand) 
            enddo        
            call random_number (xrand)        
            itur(2) = floor(xrand*num_Gauge)+1
            if(itur(2).gt.num_Gauge)itur(2) =itur(2) -1
            !put back the disturbed one
            search_list(detsel(iProc)%icho(itur(2)))%ichosen = 0
            detsel(iProc)%icho(itur(2)) = itur(1)
            call random_number (xrand)        
            itur(3) = floor(xrand*PreC)+1
            if(itur(3).gt.PreC)itur(2) =itur(3) -1
            detsel(iProc)%cita(itur(2)) = search_list(itur(1))%cita(itur(3))
            detsel(iProc)%F_matrix(itur(2),:) = search_list(itur(1))%eps(itur(3),:)            
            sM_ori = matmul(transpose(detsel(iProc)%F_matrix),detsel(iProc)%F_matrix)        
            call bsdet(sM_ori,ISTAT,det)            
            if(det.lt.singular_cri)then
                det = 0.0
            else       
                !get inverse of sM_ori, sM_ori = sM_inv
                call bssgj(sM_ori,ISTAT)
                detsel(iProc)%sM_inv = sM_ori
            endif
            detsel(iProc)%detM = det
        endif
        end subroutine disturb
        
        !singular critic conditions
        subroutine singular_critic(iProc)
        implicit double precision(a-h,o-z)
        integer itot,idif,m_sita(2)
        integer,allocatable:: icho(:)
        double precision step1,step2,step3,det,comp(2)
        double precision,allocatable:: se_vec(:,:)
        double precision,allocatable:: sM_ori(:,:)
        integer ::ISTAT,num_Gauge,PreC
        common ISTAT,num_Gauge,PreC
        tr_rad = PI/180.0
        tr_sita = 180.0/PI
        singular_cri = 0.0
        search_list(:)%ichosen = 0
        itot = size(ele,dim=2)
        allocate(icho(num_Gauge))
        allocate(se_vec(2,ISTAT))
        allocate(sM_ori(ISTAT,ISTAT))
        call random_seed () 
        do ini = 1,20
            do i=1,num_Gauge            
                call random_number (xrand)
                num = floor(xrand*itot)+1
                if(num.gt.itot)num=num-1
                if(i.lt.2)then
                    icho(i) = num
                else
                    isig = 0
                    do while(isig.lt.1)
                        do j=1,i-1
                            idif = abs(num - icho(j))
                            if(idif.lt.1)isig = 1
                        enddo
                        if(isig.lt.1)then
                            icho(i) = num
                            isig = 1
                        else
                            call random_number (xrand)
                            num = floor(xrand*itot)+1
                            if(num.gt.itot)num=num-1
                            isig = 0
                        endif
                    enddo
                endif
                search_list(icho(i))%ichosen = 1
            enddo
            do i=1,num_Gauge
                ip = icho(i)
                 tan_in = ele(1,ip)%strain(3)/(ele(1,ip)%strain(1)-ele(1,ip)%strain(2))
                 msita = 0.5*tr_sita*atan(tan_in)
                 if(msita.lt.0)then
                     m_sita(1) = msita+90
                     m_sita(2) = msita+180
                 else
                     m_sita(1) = msita
                     m_sita(2) = msita+90
                 endif
                 comp(:)=0.0
                 do j=1,2 
                     t_cita = m_sita(j)*tr_rad
                     do k = 1,ISTAT
                        step1 = 0.5*(ele(k,ip)%strain(1)+ele(k,ip)%strain(2))
                        step2 = 0.5*(ele(k,ip)%strain(1)-ele(k,ip)%strain(2))*cos(2.0*t_cita)
                        step3 = 0.5*ele(k,ip)%strain(3)*sin(2.0*t_cita)
                        se_vec(j,k) = step1+step2+step3
                        comp(j)=comp(j)+abs(se_vec(j,k))
                     enddo
                 enddo
                 if(comp(1).gt.comp(2))then
                     ic = 1
                 else
                     ic = 2
                 endif 
                 detsel(iProc)%icho(i) =  icho(i)
                 detsel(iProc)%cita(i) = m_sita(ic)
                 detsel(iProc)%F_matrix(i,:) = se_vec(ic,:)
            enddo        
            sM_ori = matmul(transpose(detsel(iProc)%F_matrix),detsel(iProc)%F_matrix)        
            call bsdet(sM_ori,ISTAT,det)
            if(det.lt.0.0)det = 0.0
            singular_cri = singular_cri+det
        enddo
        singular_cri = 0.0001*singular_cri/20.0
        end subroutine singular_critic 
        
        
        !detmax主程序，D-optimal
        subroutine detmax(iProc)
        implicit double precision(a-h,o-z)
        integer imax(2),imin
        integer,allocatable:: icho(:),cita(:)
        double precision det_max(2),det_min,sdet_t(1,1),sddd
        double precision,allocatable:: sM_inv_plus(:,:),sM_inv_y(:,:),F_matrix_plus(:,:), sM_ori(:,:)
        double precision,allocatable:: vector_y(:,:),vector_mid(:,:),vector_max(:,:),vector_min(:,:)
        integer ::ISTAT,num_Gauge,PreC
        common ISTAT,num_Gauge,PreC
        allocate(icho(num_Gauge+1),cita(num_Gauge+1))
        allocate(vector_y(ISTAT,1),vector_mid(1,ISTAT),vector_max(ISTAT,1),vector_min(ISTAT,1))
        allocate(sM_inv_plus(ISTAT,ISTAT),sM_inv_y(ISTAT,1),F_matrix_plus(num_Gauge+1,ISTAT))
        allocate(sM_ori(ISTAT,ISTAT))        
        itot = size(ele,dim=2)
        !TestTestTestTestTestTestTestTestTestTestTestTestTestTest
        !!$OMP CRITICAL
        !    print*,PreC,iTerate
        !!$OMP END CRITICAL
        !It may be noted that when attempting to estimate k loads (or MPF), 
        !the number of gauges (m) used must satisfy the relation m>k.
        !random initialization 
        findmax(iProc)%detM = 0.0
        call singular_critic(iProc)        
        call rand_init(iProc)
        do while(detsel(iProc)%detM.eq.0.0)
            call rand_init(iProc)
        enddo
        !iteration        
        iflag = 0
        iter = 1
        it_loop = 1
        do while(iflag.ne.1)            
            !det_max(1) i iteration max; det_max(2) j iteration max
            !imax(1) i iteration max; imax(2) j iteration max
            det_max(:) = 0.0
            imax(:) = 0
            !Plus one row to find detmax
            icho(1:num_Gauge) = detsel(iProc)%icho
            cita(1:num_Gauge) = detsel(iProc)%cita
            do i = 1,itot
                if(search_list(i)%ichosen.ne.0)cycle
                !isig = 0
                !loop1:do j=1,num_Gauge
                !    idif = abs(detsel(iProc)%icho(j)-i)
                !    if(idif.lt.1)then
                !        isig =1
                !        exit loop1
                !    endif
                !end do loop1
                !if(isig .gt.0)cycle
                !find detmax in one element
                do j=1,PreC
                    vector_y(:,1) = search_list(i)%eps(j,:)
                    vector_mid = matmul(transpose(vector_y),detsel(iProc)%sM_inv)
                    sdet_t = matmul(vector_mid,vector_y)
                    det_temp = sdet_t(1,1)
                    if(det_max(1).lt.det_temp)then
                        det_max(2) = det_temp
                        imax(2) = j
                        vector_max(:,1) = vector_y(:,1)
                    endif
                enddo
                !find detmax in all element
                if(det_max(1).lt.det_max(2))then
                    det_max(1) = det_max(2)
                    imax(1) = i
                endif
            enddo
            !decision of converge
            if(imax(1).eq.0)goto 99
            if(it_loop.gt.iTerate)then
                iflag = 1
                !$OMP CRITICAL
                write(10,*)">>Converged."
                    !E11.4
                write(10,"(A10,10(I10,5X))")"Gauge:",(search_list(findmax(iProc)%icho(i))%id,i=1,num_Gauge)
                write(10,"(A10,10(I10,5X))")"Theta:",(findmax(iProc)%cita(i),i=1,num_Gauge)
                write(10,*)"---------------------------------------------"
                do i=1,num_Gauge
                     write(10,"(4(ES11.4,5X))")(findmax(iProc)%F_matrix(i,j),j=1,ISTAT)
                enddo
                !$OMP END CRITICAL
            else
                !sign the search_list ichosen label to died
                search_list(imax(1))%ichosen = 1
                icho(num_Gauge+1) = imax(1)
                cita(num_Gauge+1) = search_list(imax(1))%cita(imax(2))
                !manage epslon matrix
                F_matrix_plus(1:num_Gauge,:) = detsel(iProc)%F_matrix
                F_matrix_plus(num_Gauge+1,:) =  vector_max(:,1)
                !minus  one row to find detmax_min
                sM_inv_y = matmul(detsel(iProc)%sM_inv,vector_max)
                sM_inv_plus=detsel(iProc)%sM_inv-matmul(sM_inv_y,transpose(sM_inv_y))/(1.0+det_max(1))
                det_min = 0.0
                imin = 0
                do i =1,num_Gauge+1
                    vector_y(:,1) = F_matrix_plus(i,:)
                    vector_mid = matmul(transpose(vector_y),sM_inv_plus)
                    sdet_t = matmul(vector_mid,vector_y)
                    det_temp = sdet_t(1,1)
                    if(det_min.gt.det_temp.or.i.eq.1)then
                        det_min = det_temp
                        imin = i
                        vector_min(:,1) = vector_y(:,1)
                    endif
                enddo
                isemax = icho(imin)
                !get the ordinary sM_inv and pass it to detsel
                detsel(iProc)%detM=detsel(iProc)%detM*(1.0+det_max(1))*(1.0-det_min)
                sM_inv_y = matmul(sM_inv_plus,vector_min)
                detsel(iProc)%sM_inv=sM_inv_plus+matmul(sM_inv_y,transpose(sM_inv_y))/(1.0-det_min)
                !reorganize detsel matrix
                jc=1
                do i=1,num_Gauge+1
                    if(i.ne.imin)then
                        detsel(iProc)%icho(jc) = icho(i)
                        detsel(iProc)%cita(jc) = cita(i)
                        detsel(iProc)%F_matrix(jc,:) = F_matrix_plus(i,:)
                        jc=jc+1
                    else
                        !sign the search_list ichosen label to start
                        !test find：restart makes:1\slow down the converge 2\trapped in local optimun
                        search_list(icho(i))%ichosen = 0
                        jc=jc
                    endif
                enddo                 
            endif
            !Add disturbance term to prevent local optimization
            if(abs(imax(1)-isemax).eq.0)then                
                if(findmax(iProc)%detM.lt.detsel(iProc)%detM)then
                    findmax(iProc) = detsel(iProc)
                    sM_ori = matmul(transpose(findmax(iProc)%F_matrix),findmax(iProc)%F_matrix)        
                    call bsdet(sM_ori,ISTAT,det)
                    !$OMP CRITICAL
                    write(10,*)det
                    !$OMP END CRITICAL
                endif
                it_loop = it_loop+1
99           call disturb(iProc)
                do while(detsel(iProc)%detM.eq.0.0)
                    call disturb(iProc)
                enddo
            endif
            sM_ori = matmul(transpose(detsel(iProc)%F_matrix),detsel(iProc)%F_matrix)        
            call bsdet(sM_ori,ISTAT,det)
            !$OMP CRITICAL
            if(detsel(iProc)%detM/det.gt.5)then                 
                print*,"GT than 5"
                do i=1,num_Gauge
                     write(10,"(4(ES11.4,5X))")(detsel(iProc)%F_matrix(i,j),j=1,ISTAT)
                enddo
                pause
            endif
            write(*,"(A71)")" ======================================================================"
            write(*,"(A17,T22,A3,I7,A30,I7)")"PROCESSOR","=",iProc,"OUTER LOOP ITERATION = ",iter
            write(*,"(A71)")" ----------------------------------------------------------------------"
            write(*,"(A19,T22,A3,ES11.4)")"DETERMINANT","=",det
            write(*,"(A71)")" ----------------------------------------------------------------------"
            write(*,"(A2,T12,A12,T34,A1,T50,A5,T71,A1)")"|","Guage Number","|","Angle","|"
            write(*,"(A71)")" ----------------------------------------------------------------------"
            do i=1,num_Gauge
                write(*,"(A2,T13,I7,T34,A1,T47,I7,T71,A1)")" |",search_list(detsel(iProc)%icho(i))%id,"|",detsel(iProc)%cita(i),"|"
            enddo
            write(*,"(A71)")" +-------------------------------+------------------------------------+//"
             !$OMP END CRITICAL
            !if(iter.gt.7750)pause
            iter = iter+1
        enddo
        deallocate(icho,cita,vector_y,vector_mid,vector_max,vector_min)
        deallocate(sM_inv_plus,sM_inv_y,F_matrix_plus,sM_ori)
        end subroutine detmax
        
    end module algo