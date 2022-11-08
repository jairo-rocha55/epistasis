#!/usr/bin/perl

use Parallel::ForkManager;

# inicializa el paralelizador a 8 procesos
my $pm = new Parallel::ForkManager(16);

#  vcfList is the list of files in COAD/vcf/ to process
 
open(LIST,"../COAD/vcfList");
while (<LIST>) {
  chop();
        # inicializa cada proceso a paralelizar
  $pm -> start and next;
  system("./annovarCOAD.sh $_");
  $pm -> finish;
}

$pm -> wait_all_children; 

close(LIST);
