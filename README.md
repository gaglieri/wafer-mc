;;; Sensor yield estimates by Monte Carlo

;; This program executes Monte Carlo experiments to model, simulate and analyse the
;; distributions of good and faulty tiles counts for the wafers of a sensor.
;;
;; Segments yield maps are represented as list of integers (0 or 1) representing
;; faulty or good tiles.
;; Wafers yield maps and half-layers yield maps are represented as lists of segments
;; yields maps.
;; The random generation of a segment yield map is based on a set of bernoulli
;; sample extractions with a given tile-fault-probability.
;; The counting functions receive a wafer yield map and are specialised to
;; count how many dies of a given half layer type (represented by the sublits of
;; consecutive segments yield maps) can be diced out with a fraction of faulty tiles
;; less than a given threshold.
