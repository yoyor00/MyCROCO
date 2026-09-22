#include "cppdefs.h"
#ifdef KXX 

        module read_kagrif_data
        implicit none

!$AGRIF_DO_NOT_TREAT
        character(len=*), parameter :: target_string = "KAGRIF" ! La chaîne cible à rechercher
        character(len=256) :: filename ! Nom du fichier
        character(len=1024) :: line ! Variable pour stocker chaque ligne du fichier
        character(len=1024) :: data_line ! Variable pour stocker la ligne de données
        character(len=128), dimension(:), allocatable :: data_array ! Tableau pour stocker les éléments de la ligne de données
        integer :: i, j, num_elements ! Compteurs et nombre d'éléments
        integer :: file_unit = 10 ! Unité de fichier pour la lecture
        logical :: found_target = .false. ! Indicateur pour savoir si la chaîne cible a été trouvée
        integer :: string_length
        integer :: output_unit=6
        integer :: ierr
        character(len=1) :: char
!$AGRIF_END_DO_NOT_TREAT


        contains
        
        subroutine testmod
        implicit none
        ! Vérifier si le nom du fichier est fourni en argument
c         if (command_argument_count() < 1) then
c      &    write(output_unit, *) "Usage: ./a.out <nom_du_fichier>"
c             stop 1
c         end if

        ! Récupérer le nom du fichier depuis les arguments de la ligne de commande
        call get_command_argument(1, value=filename)

        ! Ouvrir le fichier en mode lecture
        open(unit=file_unit, file=trim(filename), action="read",
     &    iostat=ierr)
        if (ierr /= 0) then
        write(output_unit, *)
     &  "Erreur lors de l ouverture du fichier : ",trim(filename)
        stop 2
        end if

    ! Boucle de lecture du fichier ligne par ligne
        do
            read(file_unit, fmt="(a)", iostat=ierr) line
            if (ierr /= 0) exit ! Sortir de la boucle si fin de fichier ou erreur de lecture

        !    Rechercher la chaîne cible
            if (index(line, target_string) > 0) then
                found_target = .true.
                ! Lire la ligne suivante
                read(file_unit, fmt="(a)", iostat=ierr) data_line
                if (ierr /= 0) then
                    write(output_unit, *)
     & "Erreur : Fin de fichier inattendue après avoir trouvé KAGRIF."
                    close(file_unit)
                    stop 3
                end if
                exit ! Sortir de la boucle après avoir trouvé et lu la ligne suivante
            end if
        end do

        close(file_unit) ! Fermer le fichier

        ! Vérifier si la chaîne cible a été trouvée
        if (.not. found_target) then
            write(output_unit, *)
     &      "La chaîne ", target_string,
     &      " n a pas été trouvée dans le fichier."
            stop 4
        end if
    
        ! Traiter la ligne de données (celle après la ligne contenant KAGRIF)
        if (found_target) then
        ! Ignorer le '#' initial et supprimer les espaces redondants en début de ligne
            data_line = trim(adjustl(replace_sharp(data_line)))
    
            ! Compter le nombre d'éléments séparés par des espaces
            num_elements = count_elements(data_line)
    
            ! Allouer de la mémoire pour le tableau de chaînes de caractères
            allocate(data_array(num_elements))
    
            ! Extraire les éléments de la ligne et les stocker dans le tableau
            call extract_elements(data_line, data_array)
    
            ! Afficher les éléments du tableau
            write(output_unit, *) "Éléments trouvés :"
            do i = 1, num_elements
                write(output_unit, *) trim(data_array(i))
            end do
    
            ! Libérer la mémoire allouée pour le tableau
            deallocate(data_array)
            end if
        
        end subroutine testmod

        ! Fonction pour compter le nombre d'éléments séparés par des espaces
        function count_elements(line) result(count)
        implicit none
        character(len=*), intent(in) :: line
        integer :: count
        integer :: i
        logical :: in_element

        count = 0
        in_element = .false.
        do i = 1, len(line)
            if (line(i:i) /= ' ' .and. .not. in_element) then
                count = count + 1
                in_element = .true.
            else if (line(i:i) == ' ') then
                in_element = .false.
            end if
        end do
        end function

        ! Subroutine pour extraire les éléments d'une ligne et les stocker dans un tableau
        subroutine extract_elements(line, array)
        implicit none
        character(len=*), intent(in) :: line
        character(len=128), dimension(:), intent(out) :: array
        integer :: i, start_index, end_index, element_index
        logical :: in_element

        start_index = 1
        element_index = 1
        in_element = .false.

        do i = 1, len(line)
            if (line(i:i) /= ' ' .and. .not. in_element) then
                start_index = i
                in_element = .true.
            else if (line(i:i) == ' ' .and. in_element) then
                end_index = i - 1
                array(element_index) = line(start_index:end_index)
                element_index = element_index + 1
                in_element = .false.
            end if
        end do

        ! Gérer le dernier élément s'il n'est pas suivi d'un espace
        if (in_element) then
            array(element_index) = line(start_index:len(line))
        end if
        end subroutine extract_elements
    
        ! Function to replace the first '#' character with a space.
        function replace_sharp(input_string) result(output_string)
        implicit none
        character(len=*), intent(in) :: input_string
        character(len=len(input_string)) :: output_string
        integer :: i
        logical :: replaced

        output_string = input_string
        replaced = .false.
        do i = 1, len(input_string)
            if (input_string(i:i) == '#' .and. .not. replaced) then
                output_string(i:i) = ' '
                replaced = .true.
            end if
        end do
        end function

        end module read_kagrif_data
!---------------------------------------------------------------------

        subroutine getlistekernels(filein,tabbkernels,nbkernels)
        use read_kagrif_data
        implicit none
        character(len=*),intent(in) :: filein
        character(len=10),dimension(1000),intent(out) :: tabbkernels
        integer,intent(out) :: nbkernels

        ! Ouvrir le fichier en mode lecture
        open(unit=file_unit, file=trim(filein), action="read",
     &    iostat=ierr)
        if (ierr /= 0) then
        write(output_unit, *) "Erreur lors de l ouverture du fichier:"
     &  , trim(filein)
        stop 2
        end if

    ! Boucle de lecture du fichier ligne par ligne
        do
            read(file_unit, fmt="(a)", iostat=ierr) line
            if (ierr /= 0) exit ! Sortir de la boucle si fin de fichier ou erreur de lecture

        !    Rechercher la chaîne cible
            if (index(line, target_string) > 0) then
                found_target = .true.
                ! Lire la ligne suivante
                read(file_unit, fmt="(a)", iostat=ierr) data_line
                if (ierr /= 0) then
                    write(output_unit, *)
     & "Erreur : Fin de fichier inattendue après avoir trouvé KAGRIF."
                    close(file_unit)
                    stop 3
                end if
                exit ! Sortir de la boucle après avoir trouvé et lu la ligne suivante
            end if
        end do

        close(file_unit) ! Fermer le fichier

        ! Vérifier si la chaîne cible a été trouvée
        if (.not. found_target) then
            write(output_unit, *)
     &      "La chaîne ", target_string,
     &      " n a pas été trouvée dans le fichier."
            stop 4
        end if
    
        ! Traiter la ligne de données (celle après la ligne contenant KAGRIF)
        if (found_target) then
        ! Ignorer le '#' initial et supprimer les espaces redondants en début de ligne
            data_line = trim(adjustl(replace_sharp(data_line)))
    
            ! Compter le nombre d'éléments séparés par des espaces
            num_elements = count_elements(data_line)
    
            ! Allouer de la mémoire pour le tableau de chaînes de caractères
            allocate(data_array(num_elements))
    
            ! Extraire les éléments de la ligne et les stocker dans le tableau
            call extract_elements(data_line, data_array)
    
            ! Afficher les éléments du tableau
            !write(output_unit, *) "Éléments trouvés :"
            nbkernels = num_elements            
            do i = 2, num_elements
                tabbkernels(i-1) = trim(data_array(i))
                !write(output_unit, *) trim(data_array(i))
            end do
            
            ! Libérer la mémoire allouée pour le tableau
            deallocate(data_array)
            end if
        end subroutine getlistekernels
              
#else
        subroutine  testmod
        end subroutine testmod
#endif        
        
#ifdef TESTROUTINE
! ifort -DKXX -DTESTROUTINE read_kxxagrif_in.F
! ./a.out AGRIF_FixedGrids.in
        program test
#ifdef KXX        
        use read_kagrif_data
        implicit none
        character(len=10),dimension( 1000 ) :: tabkernels
        integer :: nbkernels
        
        call getlistekernels('./AGRIF_FixedGrids.in',
     &      tabkernels,nbkernels)
        print *,"kernels(/zoom):",tabkernels(1:nbkernels)
#endif        
c         call testmod
        end program test
#endif        
