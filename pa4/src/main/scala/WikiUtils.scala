import java.io.Serializable
import java.util.regex.Pattern

import scala.collection.mutable.ArrayBuffer

import org.apache.hadoop.conf.Configuration
import org.apache.hadoop.io.LongWritable
import org.apache.hadoop.io.Text
import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD

import com.databricks.spark.xml.XmlInputFormat

/**
 * A Wikipedia article, parsed into a title and list of links.
 */
case class Article(title: String, links: Seq[String])

/**
 * Utility class to parse Wikipedia articles from XML format.
 */
object WikiUtils {
  /**
   * Parse a given XML file into an RDD of Article objects. Supports both paths on the local
   * file system (when running in local mode) and Amazon S3 paths starting with s3n://.
   */
  def loadArticles(sc: SparkContext, path: String): RDD[Article] = {
    val conf = new Configuration()
    conf.set(XmlInputFormat.START_TAG_KEY, "<page>")
    conf.set(XmlInputFormat.END_TAG_KEY, "</page>")
    conf.set(XmlInputFormat.ENCODING_KEY, "utf-8")
    // Do not change! These keys give read-only access to the Amazon S3 buckets where the data
    // for CS 149 is stored.
    conf.set("fs.s3n.awsAccessKeyId", "AKIAISWF5IVKSSS5Z7MQ")
    conf.set("fs.s3n.awsSecretAccessKey", "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")

    val rdd = sc.newAPIHadoopFile(
      path, classOf[XmlInputFormat], classOf[LongWritable], classOf[Text], conf)
    rdd.map(pair => parseArticle(pair._2.toString))
  }

  // Regular expressions used for finding the title, text, and link characters ([[, | and ]]).
  val TITLE_REGEX = Pattern.compile("<title>(.*?)</title>")
  val TEXT_REGEX = Pattern.compile("<text.*?>(.*?)</text>", Pattern.MULTILINE | Pattern.DOTALL)
  val LINK_REGEX = Pattern.compile("\\[\\[|]]|\\|")

  /**
   * Populate an Article object with a title and links from the XML text for its <page> element.
   */
  def parseArticle(xml: String): Article = {
    // Find the title using the title regular expression.
    val titleMatcher = TITLE_REGEX.matcher(xml)
    val title = if (titleMatcher.find()) {
      titleMatcher.group(1).toLowerCase
    } else {
      "UNTITLED"
    }

    // Find the text tag using its regular expression.
    val textMatcher = TEXT_REGEX.matcher(xml)
    val text = if (textMatcher.find()) textMatcher.group(1) else ""

    // Find the links in the text; this is a bit tricky because links in Wikipedua markup can be
    // nested, e.g. in image captions such as this:
    //
    // [[File:Colorful spring garden.jpg|thumb|180px|right|[[Spring]] flowers in April.]
    //
    // To parse them properly, we keep track of our nesting level as we read through the string.

    val links = new ArrayBuffer[String]
    val linkMatcher = LINK_REGEX.matcher(text)
    var nestingLevel = 0  // How many [[ ... ]] pairs we are inside of
    var prevOpenPos = -1  // Position after the last [[ token we saw, or -1 if last token wasn't [[
    while (linkMatcher.find()) {
      val found = linkMatcher.group()
      if (found == "[[") {
        prevOpenPos = linkMatcher.end()
        nestingLevel += 1
      } else if (found == "|") {
        if (prevOpenPos != -1) {
          links += text.substring(prevOpenPos, linkMatcher.start()).toLowerCase
        }
        prevOpenPos = -1
      } else {  // found == "]]"
        if (prevOpenPos != -1) {
          links += text.substring(prevOpenPos, linkMatcher.start()).toLowerCase
        }
        prevOpenPos = -1
        nestingLevel -= 1
      }
    }

    new Article(title, links)
  }
}
