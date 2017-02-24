name := "CS149 Spark Assignment"

version := "1.0"

scalaVersion := "2.11.7"

libraryDependencies += "org.apache.spark" %% "spark-core" % "2.1.0"

libraryDependencies += "com.databricks" %% "spark-xml" % "0.4.1"

fork := true

outputStrategy := Some(StdoutOutput)
