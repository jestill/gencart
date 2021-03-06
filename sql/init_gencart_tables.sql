-------------------------------------------------------------+
-- init_genecart_tables.sql                                  |
-- SQL Code to generate the GenCart Schema in MySQL          |
-------------------------------------------------------------+
--                                                           |
--  AUTHOR: James C. Estill                                  |
-- STARTED: 10/26/2009                                       |
-- UPDATED: 10/26/2009                                       |
--                                                           |
-- NOTES:                                                    |
--  Generated by MySQL dump 10.10                            |
--   Server version	5.0.22                               |
--                                                           |
-------------------------------------------------------------+


/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

-------------------------------+
-- tbl_family                  |
-------------------------------+
-- The following enum field for superfamily follows the code of 
-- Wicker et al. 2007 Nat Rev Genet 8(12): 973-82.
-- for LTR retrotransposons

CREATE TABLE `tbl_family` (
  `family_name_id` int(11) NOT NULL auto_increment,
  `value` text,
  `superfamily` enum('RLC','RLG','RLX') default NULL,
  PRIMARY KEY  (`family_name_id`),
  KEY `family_name_id` (`family_name_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-------------------------------+
-- tbl_generic_attribute       |
-------------------------------+

CREATE TABLE `tbl_generic_attribute` (
  `feature_id` int(11) default NULL,
  `feature_attribute_id` int(11) default NULL,
  KEY `feature_id` (`feature_id`),
  KEY `feature_attribute_id` (`feature_attribute_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

---------------------------------+
-- tbl_generic_attribute_value   |
---------------------------------+


CREATE TABLE `tbl_generic_attribute_value` (
  `feature_attribute_id` int(11) NOT NULL auto_increment,
  `value` text,
  PRIMARY KEY  (`feature_attribute_id`),
  KEY `feature_attribute_id` (`feature_attribute_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-------------------------------+
-- tbl_name                    |
-------------------------------+

CREATE TABLE `tbl_name` (
  `name_id` int(11) NOT NULL auto_increment,
  `value` text,
  PRIMARY KEY  (`name_id`),
  KEY `name_id` (`name_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-------------------------------+
-- tbl_seq_feature             |
-------------------------------+

CREATE TABLE `tbl_seq_feature` (
  `feature_id` int(11) NOT NULL auto_increment,
  `seq_id` int(11) default NULL,
  `source_id` int(11) default NULL,
  `name_id` int(11) default NULL,
  `start` int(11) default NULL,
  `midpoint` int(11) default NULL,
  `end` int(11) default NULL,
  `length` int(11) default NULL,
  `score` float default NULL,
  `strand` enum('+','-','.') default NULL,
  `frame` tinyint(4) default NULL,
  `attribute` text,
  `family_id` int(11) default NULL,
  PRIMARY KEY  (`feature_id`),
  KEY `feature_id` (`feature_id`),
  KEY `seq_id` (`seq_id`),
  KEY `source_id` (`source_id`),
  KEY `name_id` (`name_id`),
  KEY `start` (`start`),
  KEY `midpoint` (`midpoint`),
  KEY `end` (`end`),
  KEY `strand` (`strand`),
  KEY `family_id` (`family_id`),
  KEY `length` (`length`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-------------------------------+
-- tbl_seq_id                  |
-------------------------------+

CREATE TABLE `tbl_seq_id` (
  `seq_id` int(11) NOT NULL auto_increment,
  `value` text,
  `length` int(11) default NULL,
  PRIMARY KEY  (`seq_id`),
  KEY `seq_id` (`seq_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-------------------------------+
-- tbl_seq_string              |
-------------------------------+

CREATE TABLE `tbl_seq_string` (
  `seq_id` int(11) default NULL,
  `seq_string` longtext,
  KEY `seq_id` (`seq_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-------------------------------+
-- tbl_source                  |
-------------------------------+

CREATE TABLE `tbl_source` (
  `source_id` int(11) NOT NULL auto_increment,
  `value` text,
  PRIMARY KEY  (`source_id`),
  KEY `source_id` (`source_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

