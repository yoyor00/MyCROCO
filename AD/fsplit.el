(require 'fortran)

(defun main ()
  (let ((ffile (car command-line-args-left))
        (nfile (cadr command-line-args-left)))
    (with-current-buffer (find-file-noselect ffile)
      (save-excursion
        (goto-char (point-min))
        (while (not (eq (fortran-next-statement) 'last-statement))
          (fortran-fill-statement))
        (fortran-fill-statement))
      (write-file nfile))))

(when (member "-scriptload" command-line-args)
  (main))
