set nocompatible
set expandtab
set shortmess+=I
set noerrorbells visualbell t_vb=
set backupdir=/tmp//
set directory=/tmp//
set undodir=/tmp//
if has('autocmd')
  autocmd GUIEnter * set visualbell t_vb=
endif 
setlocal tabstop=4 softtabstop=4 shiftwidth=4
autocmd FileType * setlocal formatoptions-=c formatoptions-=r formatoptions-=o
autocmd BufRead,BufNewFile *.yml,*.yaml setlocal tabstop=2 softtabstop=2 shiftwidth=2
syntax on
set number
syntax enable
au BufRead,BufNewFile in.*,*.in,*.lmp set filetype=lammps
