#lang racket

;;; ITS3 yield estimates by Monte Carlo

(require math/distributions)

; configuration parameters

(define tile-fault-probability 0.01136)
(define tiles-per-segment 144)
(define segments-per-wafer 5)
(define segments-per-hl2 5)
(define segments-per-hl1 4)
(define segments-per-hl0 3)



; function to generate a random
; map of functioning and faulty tiles
; in one segment from a series of
; bernoulli experiments
(define (generate-segment tile-fault-probability)
  (map inexact->exact
       (map round
            (sample (bernoulli-dist (- 1 tile-fault-probability)) tiles-per-segment))))


; function to generate the representation
; of a wafer yield map obtained repeating
; the random experiment for one segment
; for the number of segments per wafer
(define (generate-wafer-yield-map tile-fault-probability)
  (for/list ([i (in-range segments-per-wafer)])
    (generate-segment tile-fault-probability)))

; checks that a random list
; of 0s and 1s numbers contains a number of 0s less
; than a given threshold fraction
(define (acceptable-fault-fraction? yield-list fault-fraction)
  (define num-instances (length yield-list))
  (define max-faults (inexact->exact (floor (* fault-fraction num-instances))))
  (define faults-count (count zero? yield-list))
  (cond
    [(<= faults-count max-faults) #t]
    [else #f]))

; checks that a given chip yield map has an acceptable
; fraction of faulty components
; a chip yield map is a list of lists of 0s and 1s,
; each representing the yield map of one of the
; segments composing a chip
(define (acceptable-chip? chip-yield-map fault-fraction)
  (acceptable-fault-fraction? (flatten chip-yield-map) fault-fraction))

; counts how many chips have an acceptable
; fraction of faulty components
; receives a wafer-yield-map (list of segments yield maps)
; counts chips grouping segments by a given parameter
(define (count-good-chips wafer-yield-map segments-per-chip acceptable-fault-fraction)
  (define num-chips (add1 (- (length wafer-yield-map) segments-per-chip)))
  (count (lambda (x) x)
         (for/list ([i (in-range num-chips)])
           (acceptable-chip? (take (drop wafer-yield-map i) segments-per-chip) acceptable-fault-fraction))))
    

; checks a given wafer yield map
; returns the number of hl2 device
; with acceptable tile fault fraction
; under the assumptions, can return only 1 or 0
(define (count-good-l2-chips wafer-yield-map acceptable-fault-fraction)
  (count-good-chips wafer-yield-map segments-per-hl2 acceptable-fault-fraction))

; checks a given wafer yield map
; returns the number of hl1 device
; with acceptable tile fault fraction
; under the assumptions, can return only 0, 1 or 2
(define (count-good-l1-chips wafer-yield-map acceptable-fault-fraction)
  (count-good-chips wafer-yield-map segments-per-hl1 acceptable-fault-fraction))

; checks a given wafer yield map
; returns the number of hl0 device
; with acceptable tile fault fraction
; under the assumptions, can return only 0, 1, 2 or 3
(define (count-good-l0-chips wafer-yield-map acceptable-fault-fraction)
  (count-good-chips wafer-yield-map segments-per-hl0 acceptable-fault-fraction))


; run and print experiment
(define (run-experiment num-wafers tile-fault-probability acceptable-fault-fraction)
  (printf "Number of wafers: ~a, Tile fault probability: ~a, Acceptable fault fraction: ~a\n" num-wafers tile-fault-probability acceptable-fault-fraction)
  (printf "Wafer\tHL2\tHL1\tHL0\n")
  (define yield-list
    (for/list ([i (in-range num-wafers)])
      (let ([test-wafer-yield-map (generate-wafer-yield-map tile-fault-probability)])
        (list (count-good-l2-chips test-wafer-yield-map acceptable-fault-fraction)
              (count-good-l1-chips test-wafer-yield-map acceptable-fault-fraction)
              (count-good-l0-chips test-wafer-yield-map acceptable-fault-fraction)))))
  (for ([wafer-id (in-range (length yield-list))]
        [yield-data (in-list yield-list)])
    (apply (lambda (yhl2 yhl1 yhl0) (printf "~a\t~a\t~a\t~a\n" (add1 wafer-id) yhl2 yhl1 yhl0)) yield-data))
  (newline))

  

(module+ main
  (require rackunit)

  (check-equal? (generate-segment 0.0) (make-list tiles-per-segment 1))
  (check-equal? (generate-segment 1.0) (make-list tiles-per-segment 0))

  (check-equal? (generate-wafer-yield-map 0.0) (make-list segments-per-wafer(make-list tiles-per-segment 1)))
  (check-equal? (generate-wafer-yield-map 1.0) (make-list segments-per-wafer(make-list tiles-per-segment 0)))

  (check-false (acceptable-fault-fraction? '(0 0 0 0 0 0 0 0 0 0) 0.01))
  (check-false (acceptable-fault-fraction? '(0 0 0 0 0 0 0 0 0 0) 0.1))
  (check-false (acceptable-fault-fraction? '(0 0 0 0 0 0 0 0 0 0) 0.9))
  (check-false (acceptable-fault-fraction? '(0 0 0 0 0 0 0 0 0 0) 0.99))
  (check-true  (acceptable-fault-fraction? '(0 0 0 0 0 0 0 0 0 0) 1.0))
  (check-true  (acceptable-fault-fraction? '(1 1 1 1 1 1 1 1 1 1) 0.01))
  (check-true  (acceptable-fault-fraction? '(1 1 1 1 1 1 1 1 1 1) 0.1))
  (check-true  (acceptable-fault-fraction? '(1 1 1 1 1 1 1 1 1 1) 0.4))
  (check-true  (acceptable-fault-fraction? '(1 1 1 1 1 1 1 1 1 1) 0.8))
  (check-true  (acceptable-fault-fraction? '(1 1 1 1 1 1 1 1 1 1) 0.9))
  (check-true  (acceptable-fault-fraction? '(1 1 1 1 1 1 1 1 1 1) 1.0))

  (check-false (acceptable-fault-fraction? '(0 0 1 1 1 1 1 1 1 1) 0.1))
  (check-true  (acceptable-fault-fraction? '(0 1 1 1 1 1 1 1 1 1) 0.1))
  (check-false (acceptable-fault-fraction? '(0 0 0 0 0 0 0 0 0 1) 0.8))
  (check-false (acceptable-fault-fraction? '(0 0 0 0 0 0 0 0 0 1) 0.89))
  (check-true  (acceptable-fault-fraction? '(0 0 0 0 0 0 0 0 0 1) 0.9))

  (check-false (acceptable-chip? '((0 1 1 1 1 1 1 1 1 1)
                                   (0 1 1 1 1 1 1 1 1 1)
                                   (0 1 1 1 1 1 1 1 1 1))
                                 0.095))
  (check-false (acceptable-chip? '((0 1 1 1 1 1 1 1 1 1)
                                   (0 1 1 1 1 1 1 1 1 1)
                                   (0 1 1 1 1 1 1 1 1 1))
                                 2/30))
  (check-true  (acceptable-chip? '((0 1 1 1 1 1 1 1 1 1)
                                   (0 1 1 1 1 1 1 1 1 1)
                                   (0 1 1 1 1 1 1 1 1 1))
                                 3/30))
  (check-false (acceptable-chip? '((0 1 1 1 1 1 1 1 1 1)
                                   (0 1 1 1 1 1 1 1 1 1)
                                   (0 0 1 1 1 1 1 1 1 1))
                                 0.1))
  (check-false (acceptable-chip? '((0 1 1 1 1 1 1 1 1 1)
                                   (0 0 1 1 1 1 1 1 1 1)
                                   (0 0 1 1 1 1 1 1 1 1))
                                 4/30))
  (check-true  (acceptable-chip? '((0 1 1 1 1 1 1 1 1 1)
                                   (0 0 1 1 1 1 1 1 1 1)
                                   (0 0 1 1 1 1 1 1 1 1))
                                 5/30))
  (check-true  (acceptable-chip? '((0 1 1 1 1 1 1 1 1 1)
                                   (0 0 1 1 1 1 1 1 1 1)
                                   (0 0 1 1 1 1 1 1 1 1))
                                 6/30))
  

  (check-equal? (count-good-chips '((0 1 1 1 1 1 1 1 1 1)
                                    (0 1 1 1 1 1 1 1 1 1)
                                    (0 1 1 1 1 1 1 1 1 1)) 3 0.095)
                                  0)
  (check-equal? (count-good-chips '((0 1 1 1 1 1 1 1 1 1)
                                    (0 1 1 1 1 1 1 1 1 1)
                                    (0 1 1 1 1 1 1 1 1 1)) 3 3/30)
                                  1)
  (check-equal? (count-good-chips '((0 1 1 1 1 1 1 1 1 1)
                                    (0 1 1 1 1 1 1 1 1 1)
                                    (0 1 1 1 1 1 1 1 1 1)) 3 4/30)
                                  1)
  (check-equal? (count-good-chips '((0 1 1 1 1 1 1 1 1 1)
                                    (0 1 1 1 1 1 1 1 1 1)
                                    (0 1 1 1 1 1 1 1 1 1)) 2 1/20)
                                  0)
  (check-equal? (count-good-chips '((0 1 1 1 1 1 1 1 1 1)
                                    (0 1 1 1 1 1 1 1 1 1)
                                    (1 1 1 1 1 1 1 1 1 1)) 2 1/20)
                                  1)
  (check-equal? (count-good-chips '((1 1 1 1 1 1 1 1 1 1)
                                    (0 1 1 1 1 1 1 1 1 1)
                                    (1 1 1 1 1 1 1 1 1 1)) 2 1/20)
                                  2)
  (check-equal? (count-good-chips '((1 1 1 1 1 1 1 1 1 1)
                                    (0 1 1 1 1 1 1 1 1 1)
                                    (1 1 1 1 1 1 1 1 1 1)) 2 1/30)
                                  0)
  (check-equal? (count-good-chips '((0 1 1 1 1 1 1 1 1 1)
                                    (0 1 1 1 1 1 1 1 1 1)
                                    (0 1 1 1 1 1 1 1 1 1)) 1 1/10)
                                  3)
  (check-equal? (count-good-chips '((0 1 1 1 1 1 1 1 1 1)
                                    (0 0 1 1 1 1 1 1 1 1)
                                    (0 1 1 1 1 1 1 1 1 1)) 1 1/10)
                                  2)
  (check-equal? (count-good-chips '((0 0 1 1 1 1 1 1 1 1)
                                    (0 0 1 1 1 1 1 1 1 1)
                                    (0 1 1 1 1 1 1 1 1 1)) 1 1/10)
                                  1)
  (check-equal? (count-good-chips '((0 1 1 1 1 1 1 1 1 1)
                                    (0 1 1 1 1 1 1 1 1 1)
                                    (0 1 1 1 1 1 1 1 1 1)) 1 1/20)
                                  0)
  (check-equal? (count-good-chips '((0 1 1 1 1 1 1 1 1 1)
                                    (0 0 1 1 1 1 1 1 1 1)
                                    (0 1 1 1 1 1 1 1 1 1)) 1 1/20)
                                  0)
  (check-equal? (count-good-chips '((0 0 1 1 1 1 1 1 1 1)
                                    (0 0 1 1 1 1 1 1 1 1)
                                    (0 1 1 1 1 1 1 1 1 1)) 1 1/20)
                                  0)

  (check-equal? (count-good-l2-chips (generate-wafer-yield-map 0.0) 0.01) 1)
  (check-equal? (count-good-l2-chips (generate-wafer-yield-map 0.0) 0.1) 1)
  (check-equal? (count-good-l2-chips (generate-wafer-yield-map 0.0) 0.9) 1)

  (check-equal? (count-good-l2-chips (generate-wafer-yield-map 1.0) 0.01) 0)
  (check-equal? (count-good-l2-chips (generate-wafer-yield-map 1.0) 0.1) 0)
  (check-equal? (count-good-l2-chips (generate-wafer-yield-map 1.0) 0.9) 0)

  (check-equal? (count-good-l2-chips (generate-wafer-yield-map 0.01) 0.2) 1)

  (check-equal? (count-good-l2-chips '((0 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)) 1/20)
                                  0)
  (check-equal? (count-good-l2-chips '((0 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)
                                       (1 1 1 1 1 1 1 1 1 1)
                                       (1 1 1 1 1 1 1 1 1 1)
                                       (1 1 1 1 1 1 1 1 1 1)) 1/20)
                                  1)
  (check-equal? (count-good-l2-chips '((0 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)
                                       (1 1 1 1 1 1 1 1 1 1)
                                       (1 1 1 1 1 1 1 1 1 1)) 1/20)
                                  0)

  (check-equal? (count-good-l1-chips '((0 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)) 3/40)
                                  0)

  (check-equal? (count-good-l1-chips '((0 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)
                                       (1 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)) 3/40)
                                  2)
  (check-equal? (count-good-l1-chips '((0 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)
                                       (1 1 1 1 1 1 1 1 1 1)) 3/40)
                                  1)

  (check-equal? (count-good-l0-chips '((0 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)) 2/30)
                                  0)
  (check-equal? (count-good-l0-chips '((0 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)) 3/30)
                                  3)
  (check-equal? (count-good-l0-chips '((1 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)) 2/30)
                                  1)
  (check-equal? (count-good-l0-chips '((1 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)
                                       (0 1 1 1 1 1 1 1 1 1)
                                       (1 1 1 1 1 1 1 1 1 1)) 2/30)
                                  2)

  

  (run-experiment 20 0.015 0.01)
  (run-experiment 20 0.015 0.015)
  (run-experiment 20 0.015 0.02)
  (run-experiment 20 0.020 0.02)
  

  (void))

