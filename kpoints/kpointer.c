

  main ()

{
   int i;
   char lett;

   for (i=1;i<5;++i)
     printf("\n");
   printf(" ****************   K Point Generators   *****************\n");
   printf("The programs in this directory are designed to generate special\n");
   printf("k point sets for the various Bravais lattices (2d and 3d), and\n"); 
   printf("for a one dimensional chain.  The following are included in the\n");
   printf("directory currently (as of 7/88): \n");
   for (i=1;i<2;++i)
     printf("\n");
   printf("   chain.exe    --------------- one dimensional chain\n");
   printf("   square.exe   --------------- 2-d square lattice\n");
   printf("   prect.exe    --------------- 2-d primitive rectangular lattice\n");
   printf("   crect.exe    --------------- 2-d centered rectangular lattice\n");
   printf("   hex2d.exe    --------------- 2-d hexagonal lattice \n");
   printf("   sc.exe       --------------- simple cubic lattice \n");
   printf("   fcc.exe      --------------- face-centered cubic lattice\n");
   printf("   bcc.exe      --------------- body-centered cubic lattice\n");
   printf("   ptet.exe     --------------- primitive tetragonal lattice\n");
   printf("   btet.exe     --------------- body-centered tetragonal lattice\n");
   printf("   hex3d.exe    --------------- 3-d hexagonal lattice\n");
   printf("   rhom.exe     --------------- rhombohedral (trigonal) lattice\n");
   printf("   portho.exe   --------------- primitive orthorhombic lattice\n");
   printf("   bc_ortho.exe --------------- body-centered orthorhombic lattice\n");
   printf("   oblique.exe  --------------- oblique lattice in space group P1\n");
   printf("\n");
   printf("The remaining lattices will be added to the list in due course.\n");
   printf ("Hit <CR> to continue\n");
   scanf("%c",&lett);
   printf("To use the program of your choice, simply type the name of the \n");
   printf("desired program, and it will execute.\n");
   printf("\n");
   printf("In some cases, the program will warn you about your choice of\n");
   printf("direct lattice vectors.  Please heed these warnings, as an er-\n");
   printf("roneous set of k points may result otherwise.\n");
   for (i=1;i<5;++i)
     printf("\n"); 
}
