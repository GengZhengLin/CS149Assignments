import org.apache.spark.SparkConf
import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.spark.rdd.RDD

object TopLinked {
  def main(args: Array[String]) {
    if (args.length != 2) {
      System.err.println("Usage: TopLinked <master-url> <wiki-file>")
      System.exit(1)
    }

    val master = args(0)
    val file = args(1)

    val conf = new SparkConf().setMaster(master).setAppName("TopLinked")
    val sc = new SparkContext(conf)

    /** YOUR CODE HERE **/
    // Feel free to import other Spark libraries as needed.

    sc.stop()
  }
}