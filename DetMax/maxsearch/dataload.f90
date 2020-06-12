    module dataload
    implicit double precision(a-h,o-z)
    double precision, parameter::PI=3.1415926
    !integer, parameter::PreC=180
    !integer, parameter::iTerate = 20
    type info_element
        !单元数据块，单元号，组成单元的编号，单元应变，e11,e22,e12
        integer id
        integer node(3)
        double precision strain(3)
    end type
    type info_select
        !迭代组数据块，包含迭代组矩阵F,当前Det，M=F'F的逆矩阵
        integer,allocatable:: icho(:),cita(:)
        double precision detM
        double precision,allocatable:: F_matrix(:,:)
        double precision,allocatable:: sM_inv(:,:)
    end type
    type info_search_ele
        !id代表所在的element编号
        !ichosen代表当前gauge是否被选中，若选中则禁用在当前搜索中
        !应始终保持当前矩阵的ichose为1禁用状态
        integer id
        integer ichosen
        !double precision cita(PreC)
        double precision,allocatable:: cita(:)
        double precision,allocatable:: eps(:,:)
        !double precision eps(PreC,ISTAT)        
    end type
    type info_node
        integer id
        double precision pos(3)
    end type
    type(info_element),allocatable::ele(:,:)
    type(info_node),allocatable::nod(:)
    type(info_search_ele),allocatable::search_list(:)
    type(info_select),allocatable::detsel(:),findmax(:)
    type(info_select) getmax
    double precision singular_cri
    integer iTerate
    character(len=100) rptfile,nodefile,outfile
    contains
        subroutine dataframe()
        implicit double precision(a-h,o-z)  
        character(len=100) text,abc,def,ghi,opq,fhi
        integer ::ISTAT,num_Gauge,PreC
        common ISTAT,num_Gauge,PreC
        abc = rptfile
        fhi = nodefile
        inum = 0
        istatus=1
        open(unit=2,file=fhi)
        do while(istatus.gt.0.5)
            read(2,*)text
            if(ichar(text(1:1)).gt.44.5.and.ichar(text(1:1)).lt.45.5)then
                istatus = -1
                do while(istatus.lt.-0.5)
                    read(2,*)text
                    if(ichar(text(1:1)).gt.44.5.and.ichar(text(1:1)).lt.45.5)then
                        istatus = 0
                    endif
                enddo
            endif
        enddo
        do while(istatus.lt.0.5)
            read(2,*)text
            if(ichar(text(1:1)).gt.44.5.and.ichar(text(1:1)).lt.45.5)then
                istatus = 1
            else
                inum=inum+1
            endif
        enddo
        close(2)
        allocate (nod(inum))
        open(unit=2,file=fhi)
        do while(istatus.gt.0.5)
            read(2,*)text
            if(ichar(text(1:1)).gt.44.5.and.ichar(text(1:1)).lt.45.5)then
                istatus = -1
                do while(istatus.lt.-0.5)
                    read(2,*)text
                    if(ichar(text(1:1)).gt.44.5.and.ichar(text(1:1)).lt.45.5)then
                        istatus = 0
                    endif
                enddo
            endif
        enddo        
        do i=1,inum
            read(2,*)text,nod(i)%id,nod(i)%pos(1),nod(i)%pos(2),nod(i)%pos(3)
        enddo
        !filename=trim(F_directory)//trim(ctemp)
        isnum=len_trim(abc)-4
        itr=ichar(abc(isnum:isnum))-1
        open(unit=1,file=abc)
        istatus=1
        do while(istatus.gt.0.5)
            read(1,*)text
            if(ichar(text(1:1)).gt.44.5.and.ichar(text(1:1)).lt.45.5)then
                istatus=-1
                do while(istatus.lt.-0.5)
                    read(1,*)text
                    if(ichar(text(1:1)).gt.44.5.and.ichar(text(1:1)).lt.45.5)then
                        istatus=0
                    endif
                enddo
            endif
        enddo
        nodenum=0
        do while(istatus.lt.0.5)
            read(1,*)text
            if(ichar(text(1:1)).gt.44.5.and.ichar(text(1:1)).lt.45.5)then
                istatus=1
            else
                nodenum=nodenum+1
            endif
        enddo
        close(1)
        nodenum=nodenum-1
        allocate (ele(ISTAT,nodenum))
        allocate(search_list(nodenum))
        do inn =1,nodenum
            allocate(search_list(inn)%cita(PreC),search_list(inn)%eps(PreC,ISTAT))
        enddo
        opq=abc(isnum+1:isnum+4)
        !Reading Nodeloadstatus
        do i=1,ISTAT
            if(i.gt.9)then
                m=(i-10)/10+1
                n=i-10*m
                abc(isnum:isnum)=char(itr+m)
                abc(isnum+1:isnum+1)=char(itr+n)
                abc=trim(abc(1:isnum+1))//opq
            else
                abc(isnum:isnum)=char(itr+i)
            endif
            open(unit=1,file=abc)
            def="->Loading Data:"
            ghi=trim(def)//trim(abc)
            write(*,"(A100)")ghi
            do while(istatus.gt.0.5)            
                read(1,*)text
                if(ichar(text(1:1)).gt.44.5.and.ichar(text(1:1)).lt.45.5)then
                    istatus=-1
                    do while(istatus.lt.-0.5)
                        read(1,*)text
                        if(ichar(text(1:1)).gt.44.5.and.ichar(text(1:1)).lt.45.5)then
                            istatus=0
                        endif
                    enddo
                endif
            enddo
            if(i.lt.2)then
                do k=1,nodenum
                    read(1,*)text,ele(i,k)%id,text,ele(i,k)%node(1),ele(i,k)%node(2),ele(i,k)%node(3)
                enddo
            else
               do k=1,nodenum
                    read(1,*)
               enddo
            endif
            istatus=1
            do mskip = 1,4
                read(1,*)
            enddo
            do k=1,nodenum
                read(1,*)text,text,text,text,ele(i,k)%strain(1),ele(i,k)%strain(2),text,ele(i,k)%strain(3)
                ele(i,k)%strain(1)=ele(i,k)%strain(1)*10.0**(6.0)
                ele(i,k)%strain(2)=ele(i,k)%strain(2)*10.0**(6.0)
                ele(i,k)%strain(3)=ele(i,k)%strain(3)*10.0**(6.0)
            enddo
            close(1)
        enddo
        end subroutine dataframe
    end module dataload