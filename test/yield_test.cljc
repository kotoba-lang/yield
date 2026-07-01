(ns yield-test
  (:require [clojure.test :refer [deftest is testing]]
            [yield]))
(deftest namespace-loads
  (testing "the restored CLJC namespace loads"
    (is (some? yield))))
