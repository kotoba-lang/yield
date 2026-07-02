(ns yield.monte-carlo
  "Monte Carlo simulation engine with configurable parameter distributions.
  Restored from kami-yield's `monte_carlo` module (kami-engine/kami-yield/
  src/monte_carlo.rs, deleted PR #82). Uses a deterministic LCG PRNG with
  u64 wraparound arithmetic (JVM `unchecked-*` ops on `long` — same
  2's-complement bit pattern as Rust's `wrapping_mul`/`wrapping_add`;
  unsigned interpretation via `unsigned-bit-shift-right`)."
  )

(defn monte-carlo-config [num-runs seed parameters] {:num-runs num-runs :seed seed :parameters parameters})
(defn mc-parameter [name nominal distribution] {:name name :nominal nominal :distribution distribution})

;; Distribution variants
(defn gaussian [sigma] {:kind :gaussian :sigma sigma})
(defn uniform [min max] {:kind :uniform :min min :max max})
(defn log-normal [mu sigma] {:kind :log-normal :mu mu :sigma sigma})

(defn monte-carlo-result [parameter-name values mean std-dev min max yield-pass]
  {:parameter-name parameter-name :values values :mean mean :std-dev std-dev
   :min min :max max :yield-pass yield-pass})

(defn- lcg-next-f64
  "One LCG step. Returns `[value state']`, `value` uniform in `[0,1)`."
  [state]
  (let [state' (unchecked-add (unchecked-multiply state 6364136223846793005) 1442695040888963407)
        shifted (unsigned-bit-shift-right state' 11)]
    [(/ (double shifted) (double (bit-shift-left 1 53))) state']))

(defn- lcg-next-gaussian
  "One standard-normal sample via Box-Muller. Returns `[value state']`."
  [state]
  (let [[u1 state] (lcg-next-f64 state)
        u1 (max u1 1e-15)
        [u2 state] (lcg-next-f64 state)]
    [(* (Math/sqrt (* -2.0 (Math/log u1))) (Math/cos (* 2.0 Math/PI u2))) state]))

(defn- sample
  "Sample a single value from `dist` around `nominal`. Returns `[value state']`."
  [state nominal dist]
  (case (:kind dist)
    :gaussian (let [[g state] (lcg-next-gaussian state)]
                [(+ nominal (* (:sigma dist) g)) state])
    :uniform (let [[u state] (lcg-next-f64 state)]
               [(+ (:min dist) (* (- (:max dist) (:min dist)) u)) state])
    :log-normal (let [[g state] (lcg-next-gaussian state)]
                  [(Math/exp (+ (:mu dist) (* (:sigma dist) g))) state])))

(defn- compute-stats [name values spec-min spec-max]
  (let [n (double (count values))
        mean (/ (reduce + 0.0 values) n)
        variance (/ (reduce + 0.0 (map #(Math/pow (- % mean) 2) values)) (max (- n 1.0) 1.0))
        std-dev (Math/sqrt variance)
        min-v (reduce min ##Inf values)
        max-v (reduce max ##-Inf values)
        pass-count (count (filter #(and (>= % spec-min) (<= % spec-max)) values))
        yield-pass (/ (double pass-count) n)]
    (monte-carlo-result name values mean std-dev min-v max-v yield-pass)))

(defn run-monte-carlo
  "Run Monte Carlo simulation: for each of `(:num-runs config)` runs, all
  `(:parameters config)` are sampled from their distributions and passed
  to `eval-fn` (a function of a vector of sampled values -> output value,
  checked against `[spec-min spec-max]`). Returns a vector of per-
  parameter + output `MonteCarloResult`s."
  [config eval-fn spec-min spec-max]
  (let [n (:num-runs config)
        params (:parameters config)
        np (count params)
        [param-samples output-values]
        (loop [i 0 state (:seed config)
               param-samples (vec (repeat np []))
               output-values []]
          (if (>= i n)
            [param-samples output-values]
            (let [[inputs state param-samples]
                  (reduce
                   (fn [[inputs state param-samples] [pi p]]
                     (let [[v state] (sample state (:nominal p) (:distribution p))]
                       [(conj inputs v) state (update param-samples pi conj v)]))
                   [[] state param-samples]
                   (map-indexed vector params))]
              (recur (inc i) state param-samples (conj output-values (eval-fn inputs))))))
        param-results (map-indexed (fn [i p] (compute-stats (:name p) (nth param-samples i) spec-min spec-max)) params)
        output-result (compute-stats "output" output-values spec-min spec-max)]
    (vec (concat param-results [output-result]))))
